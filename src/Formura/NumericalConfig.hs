{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE TemplateHaskell    #-}
{-# LANGUAGE TypeApplications    #-}

module Formura.NumericalConfig where

import           Cases (snakify)
import           Control.Lens
import           Control.Exception
import           Data.Aeson.TH
import           Data.Bifunctor (first)
import           Data.ByteString (ByteString)
import qualified Data.ByteString as B
import           Data.Data
import           Data.Foldable (toList)
import           Data.Maybe (fromMaybe)
import           Data.Scientific
import           Data.Text.Lens (packed)
import qualified Data.Yaml as Y

import Formura.Vec

-- Config for user input
data NumericalConfig = NumericalConfig
  { _ncLengthPerNode :: Vec Scientific
  , _ncGridPerNode :: Vec Int
  , _ncGridPerBlock :: Maybe (Vec Int)
  , _ncTemporalBlockingInterval :: Maybe Int
  , _ncFilterInterval :: Maybe Int
  , _ncMPIShape :: Maybe (Vec Int)
  , _ncWithOmp :: Maybe Int
  , _ncBoundary :: Maybe (Vec String)
  , _ncReduces :: Maybe (Vec String)
  } deriving (Eq, Ord, Read, Show, Typeable, Data)

makeClassy ''NumericalConfig


$(deriveJSON (let toSnake = packed %~ snakify
              in defaultOptions { fieldLabelModifier = toSnake . drop 3
                                , constructorTagModifier = toSnake
                                , omitNothingFields = True
                                })
  ''NumericalConfig)

data ConfigException = ConfigException String
  deriving Eq

instance Show ConfigException where
  show (ConfigException err) = err

instance Exception ConfigException

decodeConfig :: ByteString -> Either ConfigException NumericalConfig
decodeConfig str = check =<< first ConfigException (Y.decodeEither str)
  where
    -- MUSTなフィールドがあるか確認
    check :: NumericalConfig -> Either ConfigException NumericalConfig
    check cfg | null (cfg ^. ncLengthPerNode) = Left $ ConfigException "length_per_node should be a nonempty list"
              | null (cfg ^. ncGridPerNode) = Left $ ConfigException "grid_per_node should be a nonempty list"
              | maybe False null (cfg ^. ncGridPerBlock) = Left $ ConfigException "grid_per_block should be a nonempty list"
              | maybe False null (cfg ^. ncMPIShape) = Left $ ConfigException "mpi_shape should be a nonempty list"
              | maybe False (< 1) (cfg ^. ncTemporalBlockingInterval) = Left $ ConfigException "temporal_blocking_interval should be a positive integer"
              | maybe False (< 1) (cfg ^. ncFilterInterval) = Left $ ConfigException "filter_interval should be a positive integer"
              | maybe False (\i -> i `notElem` [0,1,2]) (cfg ^. ncWithOmp) = Left $ ConfigException "with_omp should be 0, 1 or 2"
              | otherwise = Right cfg

readConfig :: FilePath -> IO NumericalConfig
readConfig fn = do
  con <- B.readFile fn
  case decodeConfig con of
    Left err -> throwIO err
    Right cfg -> return cfg

-- Config for code generation
data BlockingType = NoBlocking
                  | TemporalBlocking
                      { tbGridPerBlock :: [Int]
                      , tbBlockPerNode :: [Int]
                      , tbInterval :: Int
                      }
  deriving (Eq, Ord, Read, Show, Typeable, Data)

-- | Per-axis boundary condition.  The default is periodic (torus), which is
-- the only behavior Formura had before.  Non-periodic axes are anchored
-- (drift-free) with ghost cells filled locally: as ghost slabs of the input
-- buffer without temporal blocking, and as per-sub-step overwrites of the
-- block buffers under temporal blocking (see imposeBoundaries).
data BoundaryCondition = BCPeriodic
                       | BCMirror         -- ^ Neumann: ghost mirrors the interior
                       | BCFixed Double   -- ^ Dirichlet: ghost holds a constant
  deriving (Eq, Ord, Read, Show, Typeable, Data)

-- | Global reduction declared in the yaml, e.g.
--
--     reduces: [res = absmax d, tot = sum q]
--
-- After Formura_Init and after every Formura_Forward the reduction is
-- evaluated over the whole (distributed) domain, combined with
-- MPI_Allreduce, and stored in the Formura_Navi field reduce_<name>,
-- where drivers can use it for convergence tests, diagnostics, and
-- adaptive control.
data ReduceOp = RSum | RMax | RMin | RAbsMax
  deriving (Eq, Ord, Read, Show, Typeable, Data)

parseReduce :: String -> Either ConfigException (String, ReduceOp, String)
parseReduce w = case words w of
  [name, "=", op, var] -> case op of
    "sum"    -> Right (name, RSum, var)
    "max"    -> Right (name, RMax, var)
    "min"    -> Right (name, RMin, var)
    "absmax" -> Right (name, RAbsMax, var)
    _        -> Left $ ConfigException $ "unknown reduce operator: " ++ op ++ " (expected: sum | max | min | absmax)"
  _ -> Left $ ConfigException $ "invalid reduce entry: '" ++ w ++ "' (expected: <name> = <op> <variable>)"

parseBoundaryCondition :: String -> Either ConfigException BoundaryCondition
parseBoundaryCondition w = case words w of
  ["periodic"]  -> Right BCPeriodic
  ["mirror"]    -> Right BCMirror
  ["fixed", v]  -> case reads v of
                     [(x, "")] -> Right (BCFixed x)
                     _         -> Left $ ConfigException $ "invalid fixed boundary value: " ++ v
  _ -> Left $ ConfigException $ "unknown boundary condition: '" ++ w ++ "' (expected: periodic | mirror | fixed <value>)"

data InternalConfig = InternalConfig
  { _icLengthPerNode :: [Scientific]
  , _icGridPerNode :: [Int]
  , _icSpaceInterval :: [Double]
  , _icBlockingType :: BlockingType
  , _icSleeve :: Int
  , _icSleeve0 :: Maybe Int
  , _icFilterSleeve :: Maybe Int
  , _icFilterInterval :: Maybe Int
  , _icMPIShape :: Maybe [Int]
  , _icWithOmp :: Int
  , _icBoundary :: [BoundaryCondition]
  , _icReduces :: [(String, ReduceOp, String)]
  } deriving (Eq, Ord, Read, Show, Typeable, Data)

makeClassy ''InternalConfig

defaultInternalConfig :: InternalConfig
defaultInternalConfig = InternalConfig
  { _icLengthPerNode = []
  , _icGridPerNode = []
  , _icSpaceInterval = []
  , _icBlockingType = NoBlocking
  , _icSleeve = 1
  , _icSleeve0 = Nothing
  , _icFilterSleeve = Nothing
  , _icFilterInterval = Nothing
  , _icMPIShape = Nothing
  , _icWithOmp = 0
  , _icBoundary = []
  , _icReduces = []
  }

convertConfig :: Int -> Maybe Int -> Maybe Int -> NumericalConfig -> Either ConfigException InternalConfig
convertConfig s s0 sf nc = do
  bcs <- traverse parseBoundaryCondition
           (maybe (replicate dim "periodic") toList (nc ^. ncBoundary))
  rds <- traverse parseReduce (maybe [] toList (nc ^. ncReduces))
  check (ic bcs rds)
  where
    dim = length (toList (nc ^. ncGridPerNode))
    nt = fromMaybe 1 (nc ^. ncTemporalBlockingInterval)
    totalGrids = (nc ^. ncGridPerNode) + pure (2*s*nt)
    bt = case (nc ^. ncGridPerBlock) of
           Nothing -> NoBlocking
           Just gpb ->
             case (nc ^. ncTemporalBlockingInterval) of
              Nothing -> NoBlocking
              Just nt -> let bpn = liftVec2 (div) totalGrids gpb
                          in TemporalBlocking (toList gpb) (toList bpn) nt
    ms = liftVec2 (mod) totalGrids <$> (nc ^. ncGridPerBlock)
    ic bcs rds = InternalConfig
          { _icLengthPerNode = toList $ nc ^. ncLengthPerNode
          , _icGridPerNode = toList $ nc ^. ncGridPerNode
          , _icSpaceInterval = toList $ (fmap (toRealFloat @Double) $ nc ^. ncLengthPerNode) / (fmap fromIntegral $ nc ^. ncGridPerNode)
          , _icBlockingType = bt
          , _icSleeve = s
          , _icSleeve0 = s0
          , _icFilterSleeve = sf
          , _icFilterInterval = nc ^. ncFilterInterval
          , _icMPIShape = toList <$> nc ^. ncMPIShape
          , _icWithOmp = fromMaybe 0 $ nc ^. ncWithOmp
          , _icBoundary = bcs
          , _icReduces = rds
          }
    -- 値が制約を満たすか確認
    check :: InternalConfig -> Either ConfigException InternalConfig
    check cfg | any (<0) (cfg ^. icLengthPerNode) = Left $ ConfigException "the element of length_per_node should be a positive number"
              | any (<1) (cfg ^. icGridPerNode) = Left $ ConfigException "the element of grid_per_node should be a positive integer"
              | maybe False (any (<1)) (cfg ^. icMPIShape) = Left $ ConfigException "the element of mpi_shape should be a positive integer"
              | maybe False (\ft -> ft `mod` nt /= 0) (cfg ^. icFilterInterval) = Left $ ConfigException "the filter interval is a multiple of temporal blocking interval"
              | maybe False (any (/=0)) ms = Left $ ConfigException "Inconsistent config"
              | (case cfg ^. icBlockingType of
                   TemporalBlocking gpb _ nti -> any (< 2*s*nti) gpb
                   NoBlocking -> False)
                = Left $ ConfigException "grid_per_block must be at least 2*sleeve*temporal_blocking_interval in every axis (smaller blocks are swallowed by the temporal-blocking walls and silently zero out the state)"
              | length (cfg ^. icBoundary) /= length (cfg ^. icGridPerNode) = Left $ ConfigException "boundary should list one entry per axis"
              | (case cfg ^. icBlockingType of
                   TemporalBlocking _ _ nti -> or [bc == BCPeriodic && n < 2*s*nti | (bc,n) <- zip (cfg ^. icBoundary) (cfg ^. icGridPerNode)]
                   NoBlocking -> False)
                = Left $ ConfigException "grid_per_node must be at least 2*sleeve*temporal_blocking_interval on every periodic axis (the one-sided halo of a blocked step is copied from the neighbor's grid)"
              | any (/= BCPeriodic) (cfg ^. icBoundary) && maybe False (any (> 1)) (cfg ^. icMPIShape) = Left $ ConfigException "non-periodic boundaries currently require mpi_shape [1,...]"
              | otherwise = Right cfg

nbuSize :: String -> InternalConfig -> Int
nbuSize a nc = 1 -- FIX ME
-- nbuSize a nc = head $ [ read $ drop (length kwd) opt
--                       | opt <- nc ^. ncOptionStrings
--                       , kwd `isPrefixOf` opt
--                       ] ++ [1]
--   where kwd = "nbu" ++ a

