{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE RecordWildCards  #-}
{-# LANGUAGE TypeApplications #-}
{-# LANGUAGE TypeFamilies     #-}
module Formura.Generator.Templates where

import           Control.Lens hiding (op)
import           Control.Monad
import           Data.Foldable (for_, toList)
import           Data.List
import qualified Data.Map.Strict as M
import           Data.Maybe (fromJust)
import           Text.Printf

import qualified Formura.Annotation as A
import Formura.Generator.Functions
import Formura.Generator.Types
import Formura.GlobalEnvironment
import Formura.NumericalConfig
import Formura.OrthotopeMachine.Graph
import Formura.Syntax
import Formura.Vec

scaffold :: BuildM ()
scaffold = do
  addHeader "<stdio.h>"
  addHeader "<stdlib.h>"
  addHeader "<stdbool.h>"
  addHeader "<math.h>"
  withMPI $ \_ -> addHeader "<mpi.h>"
  withOMP $ addHeader "<omp.h>"

  ic <- view (omGlobalEnvironment . envNumericalConfig)
  let gridPerNode = ic ^. icGridPerNode
  defineParam "Ns" (show $ ic ^. icSleeve)
  sequence_ [defineParam ("L" ++ show @Int i) (show v) | (i,v) <- zip [1..] gridPerNode]
  withMPI $ \mpiShape -> sequence_ [defineParam ("P" ++ show @Int i) (show v) | (i,v) <- zip [1..] mpiShape]

  defGlobalData
  navi <- defNavi
  defLocalBuffs (ic ^. icBlockingType)
  defCommBuffs

  defUtilFunctions

  defFormuraSetup
  withFirstStep defFormuraFirstStep
  withFilter defFormuraFilter
  defFormuraInit navi
  defFormuraStep
  defFormuraForward
  defFormuraFinalize

defGlobalData :: BuildM ()
defGlobalData = do
  gridPerNode <- view (omGlobalEnvironment . envNumericalConfig . icGridPerNode)
  typeName <- view (omGlobalEnvironment . gridStructTypeName)
  insntaceName <- view (omGlobalEnvironment . gridStructInstanceName)
  gridStruct <- mkFields <$> view omStateSignature
  target .= gridStruct
  gridStructType <- defGlobalTypeStruct typeName gridStruct (SoA gridPerNode)
  () <$ declGlobalVariable gridStructType insntaceName Nothing
  where
    mapType :: TypeExpr -> CType
    mapType (ElemType "Rational") = CDouble
    mapType (ElemType "int") = CInt
    mapType (ElemType "float") = CFloat
    mapType (ElemType "double") = CDouble
    mapType (ElemType x) = CRawType x
    mapType (GridType _ x) = mapType x
    mapType x = error $ "Unable to translate type to C:" ++ show x

    mkFields :: M.Map IdentName TypeExpr -> [(String, CType)]
    mkFields = M.foldrWithKey (\k t acc -> (k, mapType t):acc) []

defNavi :: BuildM [(String,String)]
defNavi = do
  navi <- mkNavi
  defGlobalTypeStruct "Formura_Navi" [(n,t) | (n,t,_) <- navi] Normal
  return [(n,v) | (n,_,v) <- navi]

mkNavi :: BuildM [(String,CType,String)]
mkNavi = do
  dim <- view (omGlobalEnvironment . dimension)
  axes <- view (omGlobalEnvironment . axesNames)
  gridPerNode <- view (omGlobalEnvironment . envNumericalConfig . icGridPerNode)
  lengthPerNode <- view (omGlobalEnvironment . envNumericalConfig . icLengthPerNode)
  spaceIntervals <- view (omGlobalEnvironment . envNumericalConfig . icSpaceInterval)
  mmpiShape <- view (omGlobalEnvironment . envNumericalConfig . icMPIShape)
  fs <- case mmpiShape of
          Just mpiShape -> do
            bases <- view (omGlobalEnvironment . commBases)
            let rs = bases ++ [map negate b | b <- bases]
                ranksTable = [(formatRank r,rank2arg r) | r <- rs ]
            return $ [ ("my_rank",CInt,"rank")
                     , ("mpi_world", CRawType "MPI_Comm", "cm")]
                  ++ [("rank_"++r,CInt,"Formura_Encode_rank"++ag) | (r,ag) <- ranksTable]
                  ++ [("offset_"++a,CInt,show l++"*i"++show i) | (a,l,i) <- zip3 axes gridPerNode [1..dim]]
                  ++ [("length_"++a,CDouble,show (l*fromIntegral m)) | (a,l,m) <- zip3 axes lengthPerNode mpiShape]
                  ++ [("total_grid_"++a,CInt,show (l*m)) | (a,l,m) <- zip3 axes gridPerNode mpiShape]
          Nothing ->
            return $ [("my_rank",CInt,"0")]
                  ++ [("offset_"++a,CInt,"0") | a <- axes]
                  ++ [("length_"++a,CDouble,show l) | (a,l) <- zip axes lengthPerNode]
                  ++ [("total_grid_"++a,CInt,show l) | (a,l) <- zip axes gridPerNode]
  reduces <- view (omGlobalEnvironment . envNumericalConfig . icReduces)
  return $ [("time_step",CInt,"0")]
        ++ [("lower_"++a,CInt,"0") | a <- axes]
        ++ [("upper_"++a,CInt,show l) | (a,l) <- zip axes gridPerNode]
        ++ [("space_interval_"++a,CDouble,show l) | (a,l) <- zip axes spaceIntervals]
        ++ [("reduce_"++nm,CDouble,"0.0") | (nm,_,_) <- reduces]
        ++ fs
  where
    r2a x | x == 0 = ""
          | x == 1 = "+1"
          | x == -1 = "-1"
          | otherwise = error "Error in r2a"
    rank2arg r = "(" ++ intercalate "," ["i" ++ show i ++ r2a x | (i,x) <- zip [1..length r] r] ++ ")"

-- NoBlocking:
--   Buff
--   Rslt
-- Temporal Blocking:
--   Tmp_Floor
--   Tmp_Walls
--   Buff
--   Rslt
defLocalBuffs :: BlockingType -> BuildM ()
defLocalBuffs NoBlocking = do
  gridStruct <- use target
  s <- maximum <$> getSleeves
  gridPerNode <- view (omGlobalEnvironment . envNumericalConfig . icGridPerNode)
  buffType <- defLocalTypeStruct "Formura_Buff" gridStruct (SoA $ map (+ (2*s)) gridPerNode)
  buff <- declLocalVariable (Just "static") buffType "buff" Nothing
  globalData <- getGlobalData
  addVariable "step_input" buff
  addVariable "step_output" globalData
  withFirstStep $ \_ -> do
    addVariable "first_input" buff
    addVariable "first_output" globalData
  withFilter $ \_ -> do
    addVariable "filter_input" buff
    addVariable "filter_output" globalData
defLocalBuffs (TemporalBlocking gridPerBlock blockPerNode nt) = do
  s <- view (omGlobalEnvironment . envNumericalConfig . icSleeve)
  dim <- view (omGlobalEnvironment . dimension)
  gridStruct <- use target
  buffType <- defLocalTypeStruct "Formura_Buff" gridStruct (SoA $ map (+ (2*s)) gridPerBlock)
  rsltType <- defLocalTypeStruct "Formura_Rslt" gridStruct (SoA gridPerBlock)
  buff <- declLocalVariable (Just "static") buffType "buff" Nothing
  rslt <- declLocalVariable (Just "static") rsltType "rslt" Nothing
  -- 床の準備
  ms <- maximum <$> getSleeves
  gridPerNode <- view (omGlobalEnvironment . envNumericalConfig . icGridPerNode)
  tmpFloorType <- defLocalTypeStruct "Formura_Tmp_Floor" gridStruct (SoA $ map (+ (2*ms)) gridPerNode)
  tmpFloor <- declLocalVariable (Just "static") tmpFloorType "tmp_floor" Nothing
  -- 壁の準備
  -- 3次元の場合の壁の大きさ
  -- x壁: [MY][MZ][NT][2*Ns][NY+2*Ns][NZ+2*Ns]
  -- y壁: [MX][MZ][NT][NX+2*Ns][2*Ns][NZ+2*Ns]
  -- z壁: [MX][MY][NT][NX+2*Ns][NY+2*Ns][2*Ns]
  for_ [([i == j | j <- [1..dim]],i) | i <- [1..dim]] $ \(flag,i) -> do
    let bs = [n | (n,b) <- zip blockPerNode flag, not b]
        gs = [if b then 2*s else n+2*s | (n,b) <- zip gridPerBlock flag]
    tmpWallType <- defLocalTypeStruct ("Formura_Tmp_Wall_" ++ show i) gridStruct (SoA $ bs ++ [nt] ++ gs)
    declLocalVariable (Just "static") tmpWallType ("tmp_wall_" ++ show i) Nothing
  globalData <- getGlobalData
  addVariable "step_input" buff
  addVariable "step_output" rslt
  withFirstStep $ \_ -> do
    addVariable "first_input" tmpFloor
    addVariable "first_output" globalData
  withFilter $ \_ -> do
    addVariable "filter_input" tmpFloor
    addVariable "filter_output" globalData

defCommBuffs :: BuildM ()
defCommBuffs = do
  bases <- view (omGlobalEnvironment . commBases)
  gridPerNode <- view (omGlobalEnvironment . envNumericalConfig . icGridPerNode)
  gridStruct <- use target
  commType <- defLocalTypeStruct "Formura_Comm_Buff" gridStruct Normal
  sleeves <- getSleeves
  for_ sleeves $ \s ->
    for_ bases $ \b -> do
      let r = formatRank b
          r' = formatRank $ map negate b
      declLocalVariable (Just "static") (CArray [if d == 1 then 2*s else n | (d,n) <- zip b gridPerNode] commType) ("send_buf" ++ show s ++ "_" ++ r) Nothing
      declLocalVariable (Just "static") (CArray [if d == 1 then 2*s else n | (d,n) <- zip b gridPerNode] commType) ("recv_buf" ++ show s ++ "_" ++ r') Nothing

defUtilFunctions :: BuildM ()
defUtilFunctions = do
  withMPI $ \mpiShape -> do
    dim <- view (omGlobalEnvironment . dimension)
    let encodeRankArg = [(CInt, "p" ++ show i) | i <- [1..dim]]
    let decodeRankArg = (CInt,"p"):[(CPtr CInt, "p" ++ show i) | i <- [1..dim]]
    defLocalFunction "Formura_Encode_rank" encodeRankArg CInt $ \_ -> case dim of
      1 -> do
        let m1:_ = mpiShape
        statement $ printf "return (p1+%d)%%%d" m1 m1
      2 -> do
        let m1:m2:_ = mpiShape
        statement $ printf "return ((p1+%d)%%%d + %d*((p2+%d)%%%d))" m1 m1 m1 m2 m2
      3 -> do
        let m1:m2:m3:_ = mpiShape
        statement $ printf "return ((p1+%d)%%%d + %d*((p2+%d)%%%d) + %d*((p3+%d)%%%d))" m1 m1 m1 m2 m2 (m1*m2) m3 m3
      _ -> error "No support"
    defLocalFunction "Formura_Decode_rank" decodeRankArg CVoid $ \_ -> case dim of
      1 -> do
        let m1:_ = mpiShape
        statement $ printf "*p1 = (int)p%%%d" m1
      2 -> do
        let m1:_:_ = mpiShape
        statement $ printf "*p1 = (int)p%%%d" m1
        statement $ printf "*p2 = (int)p/%d" m1
      3 -> do
        let m1:m2:_:_ = mpiShape
        statement $ printf "int p4 = (int)p%%%d" (m1*m2)
        statement $ printf "*p1 = (int)p4%%%d" m1
        statement $ printf "*p2 = (int)p4/%d" m1
        statement $ printf "*p3 = (int)p/%d" (m1*m2)
      _ -> error "No support"

  axes <- view (omGlobalEnvironment . axesNames)
  let toPosBody a =
        statement $ printf "return n.space_interval_%s*((i+n.offset_%s)%%n.total_grid_%s)" a a a
  for_ axes $ \a ->
    defGlobalFunction ("to_pos_"++a) [(CInt, "i"), (CRawType "Formura_Navi", "n")] CDouble (\_ -> toPosBody a)

defFormuraInit :: [(String,String)] -> BuildM ()
defFormuraInit navi = do
  defGlobalFunction "Formura_Init" [(CPtr CInt, "argc"),(CRawType "char ***", "argv"),(CRawType "Formura_Navi *","n")] CVoid (\_ -> initBody navi "MPI_COMM_WORLD")
  withMPI $ \_ ->
    () <$ defGlobalFunction "Formura_Custom_Init" [(CRawType "Formura_Navi *","n"),(CRawType "MPI_Comm", "comm")] CVoid (\_ -> initBody navi "comm")
  where
    initBody :: [(String,String)] -> String -> BuildM ()
    initBody fs comm = do
      dim <- view (omGlobalEnvironment . dimension)
      withMPI $ \mpiShape -> do
        when (comm == "MPI_COMM_WORLD") $ call "MPI_Init" ["argc","argv"]
        declScopedVariable Nothing (CRawType "MPI_Comm") "cm" (Just comm)
        declScopedVariable Nothing CInt "size" Nothing
        declScopedVariable Nothing CInt "rank" Nothing
        call "MPI_Comm_size" ["cm", "&size"]
        call "MPI_Comm_rank" ["cm", "&rank"]
        raw $ printf "if(size != %d) {\nfprintf(stderr,\"Do not match the number of MPI process!\");\nexit(1);\n}" (product mpiShape)
        statement $ "int " ++ intercalate "," ["i" ++ show i | i <- [1..dim]]
        call "Formura_Decode_rank" ("rank":["&i"++show i | i <- [1..dim]])
      -- Navi構造体の初期化
      for_ fs $ \(f,v) -> statement $ "n->" ++ f ++ " = " ++ v

      call "Formura_Setup" ("*n":replicate dim "0")
      withFirstStep $ \_ -> call "Formura_First_Step" ["n"]
      genReduces

defFormuraFirstStep :: MMGraph -> BuildM ()
defFormuraFirstStep g = do
  blockOffset <- getBlockOffsets
  buff <- getVariable "first_input"
  rslt <- getVariable "first_output"
  s <- fromJust <$> view (omGlobalEnvironment . envNumericalConfig . icSleeve0)
  defLocalFunction "Formura_First_Step_Kernel" ([(CPtr $ variableType buff, "buff"), (CPtr $ variableType rslt, "rslt"),(CRawType "Formura_Navi", "n")] ++ blockOffset) CVoid (mkKernel g s)
  defLocalFunction "Formura_First_Step" [(CRawType "Formura_Navi *", "n")] CVoid $ \_ ->
      noBlocking buff rslt s "Formura_First_Step_Kernel"

defFormuraFilter :: MMGraph -> BuildM ()
defFormuraFilter g = do
  blockOffset <- getBlockOffsets
  buff <- getVariable "filter_input"
  rslt <- getVariable "filter_output"
  s <- fromJust <$> view (omGlobalEnvironment . envNumericalConfig . icFilterSleeve)
  defLocalFunction "Formura_Filter_Kernel" ([(CPtr $ variableType buff, "buff"), (CPtr $ variableType rslt, "rslt"),(CRawType "Formura_Navi", "n")] ++ blockOffset) CVoid (mkKernel g s)
  defLocalFunction "Formura_Filter" [(CRawType "Formura_Navi *", "n")] CVoid $ \_ ->
      noBlocking buff rslt s "Formura_Filter_Kernel"

defFormuraSetup :: BuildM ()
defFormuraSetup = do
  blockOffset <- getBlockOffsets
  defLocalFunction "Formura_Setup" ([(CRawType "Formura_Navi", "n")] ++ blockOffset) CVoid  $ \_ -> setupBody
  where
    setupBody = do
      globalData <- getGlobalData
      mmg <- view omInitGraph
      mkKernel mmg 0 [globalData,globalData]

defFormuraStep :: BuildM ()
defFormuraStep = do
  blockOffset <- getBlockOffsets
  buffType <- variableType <$> getVariable "step_input"
  rsltType <- variableType <$> getVariable "step_output"
  defLocalFunction "Formura_Step" ([(CPtr buffType, "buff"), (CPtr rsltType, "rslt"),(CRawType "Formura_Navi", "n")] ++ blockOffset) CVoid stepBody
  where
    stepBody :: [CVariable] -> BuildM ()
    stepBody args = do
      s <- view (omGlobalEnvironment . envNumericalConfig . icSleeve)
      mmg <- view omStepGraph
      mkKernel mmg s args

defFormuraForward :: BuildM ()
defFormuraForward = do
  blockingType <- view (omGlobalEnvironment . envNumericalConfig . icBlockingType)
  defGlobalFunction "Formura_Forward" [(CRawType "Formura_Navi *", "n")] CVoid (\_ -> forwardBody blockingType)
  where
    forwardBody NoBlocking = do
      buff <- getVariable "step_input"
      rslt <- getVariable "step_output"
      s <- view (omGlobalEnvironment . envNumericalConfig . icSleeve)
      noBlocking buff rslt s "Formura_Step"
      withFilter $ \_ -> do
        filterInterval <- fromJust <$> view (omGlobalEnvironment . envNumericalConfig . icFilterInterval)
        raw $ printf "if (n->time_step %% %d == 0) {" filterInterval
        call "Formura_Filter" ["n"]
        raw $ "}"
      genReduces
      statement "n->time_step += 1"
    forwardBody (TemporalBlocking gridPerBlock blockPerNode nt) = do
      s <- view (omGlobalEnvironment . envNumericalConfig . icSleeve)
      dim <- view (omGlobalEnvironment . dimension)
      globalData <- getGlobalData
      tmpFloor <- getVariable "tmp_floor"
      copy globalData tmpFloor empty (repeat (2*s*nt))
      -- 通信
      rs <- isendrecv globalData (s*nt)
      -- 即時計算可能なブロックと通信待ちブロックの分離
      let b0 = [(0,m-1-d) | (n,m) <- zip gridPerBlock blockPerNode, let d = 2*s*nt `div` n]
          bs = [[(if i == j then m-1-d else 0, if i > j then m-1-d else m) | (n,m,i) <- zip3 gridPerBlock blockPerNode [1..dim], let d = 2*s*nt `div` n] | j <- [1..dim]]
          update = updateWithTB gridPerBlock blockPerNode nt
      update b0
      waitAndCopy rs tmpFloor (s*nt)
      mapM_ update bs
      copy tmpFloor globalData empty empty
      axes <- view (omGlobalEnvironment . axesNames)
      for_ axes $ \a -> statement $ printf "n->offset_%s = (n->offset_%s - %d + n->total_grid_%s)%%n->total_grid_%s" a a (s*nt) a a
      withFilter $ \_ -> do
        filterInterval <- fromJust <$> view (omGlobalEnvironment . envNumericalConfig . icFilterInterval)
        raw $ printf "if (n->time_step %% %d == 0) {" filterInterval
        call "Formura_Filter" ["n"]
        raw $ "}"
      genReduces
      statement $ "n->time_step += " ++ show nt

noBlocking :: CVariable -> CVariable -> Int -> String -> BuildM ()
noBlocking buff rslt s kernelName = do
  bcs <- view (omGlobalEnvironment . envNumericalConfig . icBoundary)
  if all (== BCPeriodic) bcs
    then noBlockingPeriodic buff rslt s kernelName
    else noBlockingAnchored buff rslt s kernelName bcs

noBlockingPeriodic :: CVariable -> CVariable -> Int -> String -> BuildM ()
noBlockingPeriodic buff rslt s kernelName = do
  globalData <- getGlobalData
  axes <- view (omGlobalEnvironment . axesNames)
  dim <- view (omGlobalEnvironment . dimension)
  -- 通信
  sendrecv globalData buff s
  copy globalData buff empty (repeat $ 2*s)
  -- 1ステップ更新
  call kernelName ([ref buff, ref rslt, "*n"] ++ replicate dim "0")
  for_ axes $ \a -> statement $ printf "n->offset_%s = (n->offset_%s - %d + n->total_grid_%s)%%n->total_grid_%s" a a s a a

-- | Emit the yaml-declared global reductions: a local accumulation loop over
-- the interior, an MPI_Allreduce (under MPI), and a store into the
-- Formura_Navi field reduce_<name>.  Runs after Formura_Init and after every
-- Formura_Forward, so drivers always see the reduction of the current state.
genReduces :: BuildM ()
genReduces = do
  reduces <- view (omGlobalEnvironment . envNumericalConfig . icReduces)
  unless (null reduces) $ do
    axes <- view (omGlobalEnvironment . axesNames)
    dim <- view (omGlobalEnvironment . dimension)
    globalData <- getGlobalData
    let gname = variableName globalData
        fields = getFields globalData
        ivs = ["i" ++ show i | i <- [1..dim]]
        loopsOpen = concat [printf "for (int %s = n->lower_%s; %s < n->upper_%s; %s++) {\n" iv a iv a iv | (iv,a) <- zip ivs axes]
        loopsClose = concat (replicate dim "}\n")
        cell v = gname ++ "." ++ v ++ concat ["[" ++ iv ++ "]" | iv <- ivs]
    for_ reduces $ \(nm, op, var) -> do
      unless (var `elem` fields) $
        error $ "reduce target '" ++ var ++ "' is not a state variable (have: " ++ show fields ++ ")"
      let (ini, upd, mpiOp) = case op of
            RSum    -> ("0.0",       "acc = acc + " ++ cell var,          "MPI_SUM")
            RMax    -> ("-INFINITY", "acc = fmax(acc, " ++ cell var ++ ")", "MPI_MAX")
            RMin    -> ("INFINITY",  "acc = fmin(acc, " ++ cell var ++ ")", "MPI_MIN")
            RAbsMax -> ("0.0",       "acc = fmax(acc, fabs(" ++ cell var ++ "))", "MPI_MAX")
      raw $ "{\ndouble acc = " ++ ini ++ ";\n"
         ++ loopsOpen ++ upd ++ ";\n" ++ loopsClose
      withMPI $ \_ ->
        raw $ printf "MPI_Allreduce(MPI_IN_PLACE, &acc, 1, MPI_DOUBLE, %s, n->mpi_world);\n" mpiOp
      raw $ "n->reduce_" ++ nm ++ " = acc;\n}"

-- | Non-periodic boundaries: an anchored, drift-free single-rank path.
-- The interior is placed symmetrically at +s on every axis, ghost slabs are
-- filled locally per axis (periodic wrap / mirror / constant), and the
-- kernel output lands back on the same logical cells, so offsets stay zero
-- and to_pos_* is the identity.  Corner ghosts become correct because each
-- axis pass sweeps the full extent of the previous axes' ghosts.
noBlockingAnchored :: CVariable -> CVariable -> Int -> String -> [BoundaryCondition] -> BuildM ()
noBlockingAnchored buff rslt s kernelName bcs = do
  globalData <- getGlobalData
  dim <- view (omGlobalEnvironment . dimension)
  ns <- view (omGlobalEnvironment . envNumericalConfig . icGridPerNode)
  let fullExt = [n + 2*s | n <- ns]
      slab a lo hi rhsOf = loopWith
        [ ( if i == a then "g" else "i" ++ show i
          , if i == a then lo else 0
          , if i == a then hi else fullExt !! i
          , 1 )
        | i <- [0 .. dim-1] ] $ \idx ->
          stack <>= [ mkIdent f buff idx @= rhsOf f (fromIdx idx) | f <- getFields buff ]
      atAxis a e comps = toIdx [ if i == a then e else c | (i, c) <- zip [0..] comps ]
  -- interior anchored at +s on every axis
  copy globalData buff empty (repeat s)
  -- ghost slabs
  for_ (zip [0 ..] (zip ns bcs)) $ \(a, (na, bc)) -> case bc of
    BCPeriodic -> do
      slab a 0 s          $ \f comps -> mkIdent f buff (atAxis a ("g+" ++ show na) comps)
      slab a (na+s) (na+2*s) $ \f comps -> mkIdent f buff (atAxis a ("g-" ++ show na) comps)
    BCMirror -> do
      slab a 0 s          $ \f comps -> mkIdent f buff (atAxis a (show (2*s-1) ++ "-g") comps)
      slab a (na+s) (na+2*s) $ \f comps -> mkIdent f buff (atAxis a (show (2*(na+s)-1) ++ "-g") comps)
    BCFixed v -> do
      slab a 0 s          $ \_ _ -> show v
      slab a (na+s) (na+2*s) $ \_ _ -> show v
  -- 1ステップ更新(出力は同じ論理セルに戻るため offset は更新しない)
  call kernelName ([ref buff, ref rslt, "*n"] ++ replicate dim "0")

updateWithTB :: [Int] -> [Int] -> Int -> [(Int,Int)] -> BuildM ()
updateWithTB gridPerBlock blockPerNode nt boundary = loopWith [("j" ++ show @Int i,l,u,1) | (i,(l,u)) <- zip [1..] boundary] $ \idx -> do
  dim <- view (omGlobalEnvironment . dimension)
  s <- view (omGlobalEnvironment . envNumericalConfig . icSleeve)
  tmpFloor <- getVariable "tmp_floor"
  tmpWalls <- forM [([i == j | j <- [1..dim]],i) | i <- [1..dim]] $ \(flag,i) -> do
    let gs = [if b then 2*s else n+2*s | (n,b) <- zip gridPerBlock flag]
    tmpWall <- getVariable ("tmp_wall_" ++ show i)
    return (flag, gs, tmpWall)
  buff <- getVariable "step_input"
  rslt <- getVariable "step_output"
-- - 床の読み込み
  let floorOffset = toIdx ["+" ++ show (n*(m-1)) ++ "-" ++ show n ++ "*" | (n,m) <- zip gridPerBlock blockPerNode] <> idx
  copy tmpFloor rslt floorOffset empty
-- - NT段更新
  loopWith [("it",0,nt,1)] $ \it -> do
--   - buffに床をセット
    copy rslt buff empty empty
--   - buffに壁をセット
    for_ tmpWalls $ \(flag, gs, tmpWall) -> do
      let idx0 = (toIdx [i | (i,b) <- zip (fromIdx idx) flag, not b]) >< it
      loop gs $ \idx' ->
        stack <>= [mkIdent f buff (idx' <> toIdx [if b then n else 0 | (b,n) <- zip flag gridPerBlock]) @= mkIdent f tmpWall (idx0 >< idx') | f <- (getFields tmpWall), f `elem` (getFields buff)]
      return ()
--   - 1段更新
--     The block frame: the floor holds cells shifted by 2*s*nt (the copy
--     into tmp_floor above), so buff slot b of this block holds, at
--     sub-step it, the cell floorOffset + b + it*s - 2*s*nt of the state
--     array; n->offset_* is advanced only once per Formura_Forward (by
--     s*nt, below).  The LoadIndex emission in mkKernel assumes slot b
--     holds cell b - 2*s + block_offset, so the displacement handed to the
--     kernel is floorOffset + s*(it + 2 - 2*nt): it grows by s per
--     sub-step and equals floorOffset - s*(nt - 1) at the last one.
    call "Formura_Step" ([ref buff, ref rslt,"*n"] ++ [o ++ "+" ++ show s ++ "*(it+" ++ show (2 - 2*nt) ++ ")" | o <- fromIdx floorOffset])
--   - 壁の書き出し
    for_ tmpWalls $ \(flag, gs, tmpWall) -> do
      let idx0 = (toIdx [i | (i,b) <- zip (fromIdx idx) flag, not b]) >< it
      loop gs $ \idx' ->
        stack <>= [mkIdent f tmpWall (idx0 >< idx')  @= mkIdent f buff idx' | f <- (getFields buff), f `elem` (getFields tmpWall)]
      return ()
-- - 床の書き出し
  copy rslt tmpFloor empty floorOffset

defFormuraFinalize :: BuildM ()
defFormuraFinalize =
  defGlobalFunction "Formura_Finalize" [] CVoid $ \_ ->
    withMPI $ \_ -> call "MPI_Finalize" []

data MMRange = MMRange
  { lower :: !Int
  , upper :: !Int
  }

instance Semigroup MMRange where
  x <> y = MMRange { lower = lower x `min` lower y
                   , upper = upper x `max` upper y
                   }

instance Monoid MMRange where
  mempty = MMRange 0 0

mergeRange :: MMRange -> MMRange -> MMRange
mergeRange x y = MMRange { lower = lower x + lower y
                         , upper = upper x + upper y
                         }

fromCursor :: Vec Int -> MMRange
fromCursor s = MMRange { lower = minimum s
                       , upper = maximum s
                       }

toSize :: MMRange -> Int
toSize (MMRange {..}) = abs lower + abs upper

toOffset :: MMRange -> Int
toOffset = abs . lower

calcRange :: MMGraph -> M.Map OMNodeID MMRange
calcRange = M.foldlWithKey (\acc k (Node mi _ _) -> M.insert k (worker mi acc) acc) M.empty
  where
    worker :: MMInstruction -> M.Map OMNodeID MMRange -> MMRange
    worker mi tbl = M.foldl' go mempty mi
      where go rng (Node (LoadCursorStatic s _) _ _) = rng <> fromCursor s
            go rng (Node (LoadCursor s oid) _ _) = rng <> mergeRange (fromCursor s) (tbl M.! oid)
            go rng _ = rng

mkKernel :: MMGraph -> Int -> [CVariable] -> BuildM ()
mkKernel mmg sleeve args = do
  -- FIX ME
  bcs <- view (omGlobalEnvironment . envNumericalConfig . icBoundary)
  -- An all-periodic program runs on the shifting frame: the surrounding
  -- machinery copies the state at 2*sleeve and moves n->offset_* by sleeve
  -- every step.  Any declared wall anchors the frame instead: the copy
  -- margin is sleeve and the offsets stay fixed.  The grid-index emission
  -- below must subtract whichever margin the copies used.
  let copyMargin | all (== BCPeriodic) bcs = 2*sleeve
                 | otherwise = sleeve
  let outputSize = getSize $ variableType $ args !! 1
      inputSize = map (+ (2*sleeve)) outputSize
  -- Every node that stores to an output variable must use the same range
  -- (-sleeve, +sleeve) regardless of its own stencil extent.  The kernel
  -- writes results at the raw loop index while reading inputs at
  -- idx + toOffset(range), so the range offset is the per-step shift of
  -- that variable inside the array.  The surrounding machinery (halo
  -- copies, n->offset_* updates) assumes one uniform shift of `sleeve`
  -- cells per step for all variables; per-node ranges break that contract
  -- whenever stencil radii differ between variables (e.g. a pointwise
  -- update next to a radius-1 update), silently corrupting every
  -- cross-variable read from the second step on.  Intermediate (non-void)
  -- nodes keep their own ranges: their consumers compensate via
  -- toOffset(consumer) - toOffset(producer), which stays non-negative
  -- because consumer ranges only grow.
  let globalRange = MMRange { lower = negate sleeve, upper = sleeve }
      adjustRange k r = case M.lookup k mmg of
        Just (Node _ (ElemType "void") _) -> globalRange
        _ -> r
      rangeTable = M.mapWithKey adjustRange (calcRange mmg)
  axes <- view (omGlobalEnvironment . axesNames)
  -- 中間変数を生成するかどうかを判定する
  -- 型が void なら最後に Store されているはずなので、中間変数を生成しない
  -- そうでないなら、最後に型に合う変数に結果を書き込む必要がある (最後に命令として Store がない)
  let isVoid (ElemType "void") = True
      isVoid _                 = False

      tmpPrefix = "formura_mn"
      tmpName oid = tmpPrefix ++ show oid

      insertMNStore :: OMNodeID -> MMInstruction -> MMInstruction
      insertMNStore oid mm = M.insert (MMNodeID mmid) node mm
        where mmid = M.size mm
              node = Node (Store (tmpName oid) (MMNodeID $ mmid-1)) (ElemType "void") []
  let (tmps, mmg') = M.foldrWithKey (\omid (Node mm omt x) (ts,g) -> if isVoid omt then (ts,M.insert omid (Node mm omt x) g) else ((omid,omt):ts,M.insert omid (Node (insertMNStore omid mm) omt x) g)) ([],M.empty) mmg
  -- 中間変数の宣言
  let mapTmpType :: Int -> OMNodeType -> CType
      mapTmpType _ (ElemType "Rational") = CDouble
      mapTmpType _ (ElemType "void") = CVoid
      mapTmpType _ (ElemType "double") = CDouble
      mapTmpType _ (ElemType "float") = CFloat
      mapTmpType _ (ElemType "int") = CInt
      mapTmpType _ (ElemType x) = CRawType x
      mapTmpType ds (GridType _ x) = CArray [s-ds | s <- inputSize] (mapTmpType ds x)
      mapTmpType _ _ = error "Invalid type"
  sequence_ [declScopedVariable (Just "static") (mapTmpType (toSize $ rangeTable M.! omid) omt) (tmpName omid) Nothing | (omid,omt) <- tmps]
  -- 命令の変換
  let formatNode i = "a" ++ show i
      mapType' :: MicroNodeType -> CType
      mapType' (ElemType "Rational") = CDouble
      mapType' (ElemType "void")     = CVoid
      mapType' (ElemType "double")   = CDouble
      mapType' (ElemType "float")    = CFloat
      mapType' (ElemType "int")      = CInt
      mapType' (ElemType x)          = CRawType x
      mapType' _                     = error "Invalid type"
  let genMicroInst idx _ _ (Store n x) _ _ | tmpPrefix `isPrefixOf` n = stack <>= [(n ++ show idx) @= formatNode x]
                                           | otherwise = stack <>= [(mkIdent n (args !! 1) idx) @= formatNode x]
      genMicroInst idx rng mmid mi mt annot =
        let decl x = declScopedVariable Nothing (mapType' mt) (formatNode mmid) (Just x) >> return ()
        in decl $ case mi of
            (LoadCursorStatic d n) -> mkIdent n (args !! 0) (idx <> (toIdx . toList $ d + pure (toOffset rng)))
            (LoadCursor d oid) -> tmpName oid ++ show (idx <> (toIdx . toList $ d + pure (toOffset rng - (toOffset $ rangeTable M.! oid))))
            (Imm r) -> show (realToFrac r :: Double)
            (Uniop op a) | "external-call" `isPrefixOf` op -> (fromJust $ stripPrefix "external-call/" op) ++ "(" ++ formatNode a ++ ")"
                         | otherwise -> op ++ formatNode a
            (Binop op a b) | op == "**" -> "pow(" ++ formatNode a ++ "," ++ formatNode b ++ ")"
                           | otherwise -> formatNode a ++ op ++ formatNode b
            (Triop _ a b c) -> formatNode a ++ "?" ++ formatNode b ++ ":" ++ formatNode c
            -- The grid index must follow the same displacement as the
            -- array loads beside it.  A node's reads sit at
            -- idx + d + toOffset rng in buffers whose slot b holds the
            -- cell b - copyMargin (+ n.offset on the shifting frame), so
            -- the cell this iteration evaluates is
            --   idx + cursor + toOffset rng - copyMargin + n.offset + block_offset,
            -- with the Shift cursor of an inlined instance recorded in its
            -- MMLocation annotation.  block_offset_* is the displacement of
            -- the kernel's local frame from the state array: zero on the
            -- plain path, and under temporal blocking the block's base
            -- corrected by the floor copy margin and the sub-step
            -- (updateWithTB), because the content moves by the sleeve at
            -- every sub-step while n.offset is advanced once per
            -- Formura_Forward.  Init kernels run at
            -- sleeve zero with zero-range stores, so their emission stays
            -- the historical idx + n.offset form.  The sum is reduced
            -- modulo the global grid size exactly as to_pos_* does: on the
            -- shifting frame the raw sum leaves [0, total) once n.offset
            -- has wrapped, which would hand a non-periodic coefficient the
            -- position y + L instead of y.
            LoadIndex i ->
              let cursorShift = case A.viewMaybe annot of
                    Just (MMLocation _ c) -> toList c !! i
                    Nothing -> 0 :: Int
                  correction = cursorShift + toOffset rng - copyMargin
                  axis = axes !! i
                  total = "n.total_grid_" ++ axis
                  unwrapped = (fromIdx idx) !! i
                    ++ (if correction == 0
                          then ""
                          else "+(" ++ show correction ++ ")")
                    ++ "+n.offset_" ++ axis ++ "+block_offset_" ++ show (i+1)
              in "((" ++ unwrapped ++ ")%" ++ total ++ "+" ++ total ++ ")%" ++ total
            Naryop op xs | Just f <- stripPrefix "external-call/" op ->
              f ++ "(" ++ intercalate "," (map formatNode xs) ++ ")"
            x -> error $ "Unimplemented for keyword: " ++ show x
  let genMMInst :: MMRange -> MMInstruction -> BuildM ()
      genMMInst rng mm = loop [s'- toSize rng | s' <- inputSize] $ \idx -> sequence_ [genMicroInst idx rng mmid mi mt annot | (mmid, Node mi mt annot) <- M.toAscList mm]
  sequence_ [genMMInst (rangeTable M.! omid) mm | (omid, Node mm _ _) <- M.toAscList mmg']
