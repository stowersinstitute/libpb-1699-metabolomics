-- https://stackoverflow.com/questions/17287080/multi-input-multi-output-compilers-with-shake
-- https://stackoverflow.com/questions/53566429/shake-how-to-set-an-environment-variable-upon-invocation
-- https://mmhaskell.com/blog/2017/3/6/making-sense-of-multiple-monads
-- https://github.com/ndmitchell/shake/issues/170
-- https://github.com/ndmitchell/shake/issues/430#issuecomment-196054985
-- https://stackoverflow.com/questions/35938956/haskell-shake-special-rule-for-building-directories


-- import Prelude ($)
import Development.Shake hiding ((*>))
-- import Development.Shake
import Development.Shake.FilePath
-- import qualified System.Environment as Env

main = shakeArgs shakeOptions $ do
    want [ "out/fig/pca-categorized-primary.pdf"
         , "out/fig/pca-categorized-lipids.pdf"
--     ,"out"</>"lipids"</>"opls"</>"no-outlier"</>"Liver"</>"Ref"</>"PvT.csv"
         ]
    "out/fig/pca-categorized-primary.pdf" %> \out -> do
      need ["src/python/pca-categorized-primary.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-primary.py --astyanax /run/media/nine/SamsungPortable/stowers/data/Astyanax/metabolome/jenna-metabolomics.csv --mammals-annotation /run/media/nine/SamsungPortable/stowers/data/Mammals/ma-metabolme/mammalian-metabolome-annotation.csv --mammals-normalized /run/media/nine/SamsungPortable/stowers/data/Mammals/ma-metabolme/mammalian-metabolome-normalized.csv --compounds /run/media/nine/SamsungPortable/stowers/data/databases/kegg/compounds.json --hmdb /run/media/nine/SamsungPortable/stowers/data/databases/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/pca-categorized-primary.pdf"

    "out/fig/pca-categorized-lipids.pdf" %> \out -> do
--       features <- getEnvWithDefault "" "PY38_ENS_FEATURES"
--       https://stackoverflow.com/questions/35938956/haskell-shake-special-rule-for-building-directories
--       needDir [takeDirectory out];
--       cmd "make -p" [takeDirectory out];
      need ["src/python/pca-categorized-lipids.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-lipids.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --exclude-outlier True --output ./out/fig/pca-categorized-lipids.pdf"
--     "out"</>"lipids"</>"opls"</>"no-outlier"</>"Liver"</>"Ref"</>"PvT.csv" *> \out -> do
--       cmd_ ("python3 ./src/python/orthogonal-lipids-pop-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json")
