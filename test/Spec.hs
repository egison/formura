module Main (main) where

import Test.Hspec (hspec)
import qualified NumericalConfigSpec

main :: IO ()
main = hspec NumericalConfigSpec.spec
