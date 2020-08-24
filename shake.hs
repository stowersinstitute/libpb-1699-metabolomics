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
    want [ "out/fig/no-outliers/pca-categorized-primary.pdf"
         , "out/fig/no-outliers/pca-categorized-lipids.pdf"
         , "out/fig/no-outliers/pca-categorized-legend.pdf"
         , "out/fig/no-outliers/orotic-acid-plot.pdf"
         , "out/fig/no-outliers/sugar-phosphate-plot.pdf"
         , "out/fig/no-outliers/inflammation-plot.pdf"
         , "out/work/primary/opls/no-outliers/Nucleotides/Muscle/Ref/PvT.csv"
--     ,"out"</>"lipids"</>"opls"</>"no-outlier"</>"Liver"</>"Ref"</>"PvT.csv"
         ]

    -- categorized pca for primary
    "out/fig/no-outliers/pca-categorized-primary.pdf" %> \out -> do
      need ["src/python/pca-categorized-primary.py"]
--       PIPENV_VENV_IN_PROJECT
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-primary.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds /run/media/nine/SamsungPortable/stowers/data/databases/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/pca-categorized-primary.pdf"

    -- categorized pca plot for lipids
    "out/fig/no-outliers/pca-categorized-lipids.pdf" %> \out -> do
--       features <- getEnvWithDefault "" "PY38_ENS_FEATURES"
--       https://stackoverflow.com/questions/35938956/haskell-shake-special-rule-for-building-directories
--       needDir [takeDirectory out];
--       cmd "make -p" [takeDirectory out];
      need ["src/python/pca-categorized-lipids.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-lipids.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --exclude-outlier True --output ./out/fig/no-outliers/pca-categorized-lipids.pdf"

    -- legend for categorized pca plots
    "out/fig/no-outliers/pca-categorized-legend.pdf" %> \out -> do
--       features <- getEnvWithDefault "" "PY38_ENS_FEATURES"
--       https://stackoverflow.com/questions/35938956/haskell-shake-special-rule-for-building-directories
--       needDir [takeDirectory out];
--       cmd "make -p" [takeDirectory out];
      need ["src/python/pca-categorized-legend.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-legend.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds /run/media/nine/SamsungPortable/stowers/data/databases/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/pca-categorized-legend.pdf"

    -- orotic acid plot
    "out/fig/no-outliers/orotic-acid-plot.pdf" %> \out -> do
      need ["src/python/orotic-acid-plot.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/orotic-acid-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds /run/media/nine/SamsungPortable/stowers/data/databases/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/orotic-acid-plot.pdf"

    -- sugar phosphate plot
    "out/fig/no-outliers/sugar-phosphate-plot.pdf" %> \out -> do
      need ["src/python/sugar-phosphate-plot.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/sugar-phosphate-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds /run/media/nine/SamsungPortable/stowers/data/databases/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/sugar-phosphate-plot.pdf"

    -- proinflammatory metabolite plot
    "out/fig/no-outliers/inflammation-plot.pdf" %> \out -> do
      need ["src/python/sugar-phosphate-plot.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/inflammation-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds /run/media/nine/SamsungPortable/stowers/data/databases/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/inflammation-plot.pdf"

    -- OPLS on primary data
    "out/work/primary/opls/no-outliers/Nucleotides/Muscle/Ref/PvT.csv" %> \out -> do
      need ["src/python/opls/primary-pop-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/primary-pop-compare.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet /opt/src/stowers/pipelines/jenna-metabolomics/data/primary/sample-sheet.csv --hmdb ./data/hmdb/hmdb.json --compounds ./data/kegg/compounds.json --output-dir out/work/primary"
--     "out"</>"lipids"</>"opls"</>"no-outlier"</>"Liver"</>"Ref"</>"PvT.csv" *> \out -> do
--       cmd_ ("python3 ./src/python/orthogonal-lipids-pop-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json")
