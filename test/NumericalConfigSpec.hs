{-# LANGUAGE OverloadedStrings #-}
module NumericalConfigSpec (spec) where

import Control.Lens (view)
import qualified Data.ByteString.Char8 as B
import Data.Either (isLeft)
import Test.Hspec

import Formura.Vec
import Formura.NumericalConfig

spec :: Spec
spec = do
  describe "Valid case" $ do
    it "1d config with Temporal blocking" $ do
      let cfg = B.unlines [ "length_per_node: [1.0]"
                          , "grid_per_node: [10]"
                          , "grid_per_block: [10]"
                          , "mpi_shape: [1]"
                          , "temporal_blocking_interval: 5"
                          ]
          cfg' = NumericalConfig
                  { _ncLengthPerNode = Vec [1.0]
                  , _ncGridPerNode = Vec [10]
                  , _ncMPIShape = Just $ Vec [1]
                  , _ncGridPerBlock = Just (Vec [10])
                  , _ncTemporalBlockingInterval = Just 5
                  , _ncFilterInterval = Nothing
                  , _ncWithOmp = Nothing
                  , _ncBoundary = Nothing
                  , _ncReduces = Nothing
                  }
      decodeConfig cfg `shouldBe` (Right cfg')
    it "3d config with Temporal blocking" $ do
      let cfg = B.unlines [ "length_per_node: [1.0,2.0,3.0]"
                          , "mpi_shape: [2,2,2]"
                          , "temporal_blocking_interval: 5"
                          , "grid_per_node: [10,10,10]"
                          , "grid_per_block: [10,10,10]"
                          , "filter_interval: 100"
                          , "with_omp: 1"
                          ]
          cfg' = NumericalConfig
                  { _ncLengthPerNode = Vec [1.0,2.0,3.0]
                  , _ncGridPerNode = Vec [10,10,10]
                  , _ncMPIShape = Just $ Vec [2,2,2]
                  , _ncGridPerBlock = Just (Vec [10,10,10])
                  , _ncTemporalBlockingInterval = Just 5
                  , _ncFilterInterval = Just 100
                  , _ncWithOmp = Just 1
                  , _ncBoundary = Nothing
                  , _ncReduces = Nothing
                  }
      decodeConfig cfg `shouldBe` (Right cfg')
    it "3d config without Temporal blocking" $ do
      let cfg = B.unlines [ "length_per_node: [1.0,2.0,3.0]"
                          , "mpi_shape: [2,2,2]"
                          , "grid_per_node: [10,10,10]"
                          ]
          cfg' = NumericalConfig
                  { _ncLengthPerNode = Vec [1.0,2.0,3.0]
                  , _ncGridPerNode = Vec [10,10,10]
                  , _ncMPIShape = Just $ Vec [2,2,2]
                  , _ncGridPerBlock = Nothing 
                  , _ncTemporalBlockingInterval = Nothing
                  , _ncFilterInterval = Nothing
                  , _ncWithOmp = Nothing
                  , _ncBoundary = Nothing
                  , _ncReduces = Nothing
                  }
      decodeConfig cfg `shouldBe` (Right cfg')
  describe "Blocking with walls" $ do
    let converted b = convertConfig 1 Nothing Nothing NumericalConfig
          { _ncLengthPerNode = Vec [1.0,1.0,0.4]
          , _ncGridPerNode = Vec [16,16,4]
          , _ncMPIShape = Nothing
          , _ncGridPerBlock = Just (Vec [22,22,10])
          , _ncTemporalBlockingInterval = Just 3
          , _ncFilterInterval = Nothing
          , _ncWithOmp = Nothing
          , _ncBoundary = Just (Vec b)
          , _ncReduces = Nothing
          }
    it "accepts temporal blocking with walls on a single rank" $ do
      fmap (view icBlockingType) (converted ["mirror","fixed 0.0","fixed 1.5"])
        `shouldBe` Right (TemporalBlocking [22,22,10] [1,1,1] 3)
    it "rejects a periodic axis shorter than the one-sided halo of a blocked step" $ do
      converted ["mirror","fixed 0.0","periodic"] `shouldSatisfy` isLeft
    it "accepts the same axis without temporal blocking" $ do
      let cfg = B.unlines [ "length_per_node: [1.0,1.0,0.4]"
                          , "grid_per_node: [16,16,4]"
                          , "boundary: [mirror, fixed 0.0, periodic]"
                          ]
      fmap (view icBlockingType) (convertConfig 1 Nothing Nothing =<< decodeConfig cfg)
        `shouldBe` Right NoBlocking
  describe "Invalid case" $ do
    it "don't exist MUST fields" $ do
      let cfg = B.unlines [ "length_per_node: [1.0]"
                          , "grid_per_block: [10]"
                          ]
      decodeConfig cfg `shouldSatisfy` isLeft
    it "has empty list" $ do
      let cfg = B.unlines [ "length_per_node: [1.0]"
                          , "grid_per_node: [10]"
                          , "grid_per_block: []"
                          , "mpi_shape: []"
                          , "temporal_blocking_interval: 5"
                          ]
      decodeConfig cfg `shouldSatisfy` isLeft
    it "has a negative temporal_blocking_interval" $ do
      let cfg = B.unlines [ "length_per_node: [1.0]"
                          , "grid_per_node: [10]"
                          , "grid_per_block: [10]"
                          , "mpi_shape: [1]"
                          , "temporal_blocking_interval: -5"
                          ]
      decodeConfig cfg `shouldSatisfy` isLeft
