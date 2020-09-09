-- https://stackoverflow.com/questions/17287080/multi-input-multi-output-compilers-with-shake
-- https://stackoverflow.com/questions/53566429/shake-how-to-set-an-environment-variable-upon-invocation
-- https://mmhaskell.com/blog/2017/3/6/making-sense-of-multiple-monads
-- https://github.com/ndmitchell/shake/issues/170
-- https://github.com/ndmitchell/shake/issues/430#issuecomment-196054985
-- https://stackoverflow.com/questions/35938956/haskell-shake-special-rule-for-building-directories


import Development.Shake hiding ((*>))
import Development.Shake.FilePath

main = shakeArgs shakeOptions $ do
    want [ "out/fig/no-outliers/pca-categorized-primary.pdf"
         , "out/fig/no-outliers/pca-categorized-lipids.pdf"
         , "out/fig/no-outliers/pca-categorized-legend.pdf"
         , "out/fig/no-outliers/orotic-acid-plot.pdf"
         , "out/fig/no-outliers/sugar-phosphate-plot.pdf"
         , "out/fig/no-outliers/metabolites-of-interest.pdf"
         , "out/work/primary/opls/no-outliers/Nucleotides/Muscle/Ref/PvT.csv"
         , "out/work/primary/glm/singlefactor/no-outliers/Nucleotides/Muscle/Ref/PvT.csv"
         , "out/supp/no-outliers/primary-pop-compare-significance.xlsx"
         , "out/supp/outliers/primary-pop-compare-significance.xlsx"
         , "out/fig/no-outliers/ratios-combined.pdf"
         , "out/work/lipids/opls/outliers/category/Sphingolipids/Muscle/Ref/PvT.csv" -- Lipid OPLS
         , "out/work/lipids/glm/singlefactor/no-outliers/Muscle/Ref/PvT.csv" -- Lipid GLM
         , "out/work/lipidcats/opls/Liver/Ref/Classes/PvT.csv" -- Lipid cats/classes OPLS
         ]

    -- categorized pca for primary
    "out/fig/no-outliers/pca-categorized-primary.pdf" %> \out -> do
      need ["src/python/pca-categorized-primary.py"]
--       PIPENV_VENV_IN_PROJECT
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-primary.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/pca-categorized-primary.pdf"

    -- categorized pca plot for lipids
    "out/fig/no-outliers/pca-categorized-lipids.pdf" %> \out -> do
      need ["src/python/pca-categorized-lipids.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-lipids.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --exclude-outlier True --output ./out/fig/no-outliers/pca-categorized-lipids.pdf"

    -- legend for categorized pca plots
    "out/fig/no-outliers/pca-categorized-legend.pdf" %> \out -> do
      need ["src/python/pca-categorized-legend.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-legend.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/pca-categorized-legend.pdf"

    -- orotic acid plot
    "out/fig/no-outliers/orotic-acid-plot.pdf" %> \out -> do
      need ["src/python/orotic-acid-plot.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/orotic-acid-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/orotic-acid-plot.pdf"

    -- sugar phosphate plot
    "out/fig/no-outliers/sugar-phosphate-plot.pdf" %> \out -> do
      need ["src/python/sugar-phosphate-plot.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/sugar-phosphate-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet /opt/src/stowers/pipelines/jenna-metabolomics/data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/sugar-phosphate-plot.pdf"

    -- proinflammatory metabolite plot
    "out/fig/no-outliers/metabolites-of-interest.pdf" %> \out -> do
      need ["src/python/metabolites-of-interest.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/metabolites-of-interest.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet /opt/src/stowers/pipelines/jenna-metabolomics/data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier True --output ./out/fig/no-outliers/metabolites-of-interest.pdf"

    -- OPLS primary data
    "out/work/primary/opls/no-outliers/Nucleotides/Muscle/Ref/PvT.csv" %> \out -> do
      need ["src/python/opls/primary-pop-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/primary-pop-compare.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet /opt/src/stowers/pipelines/jenna-metabolomics/data/primary/sample-sheet.csv --hmdb ./data/hmdb/hmdb.json --compounds ./data/kegg/compounds.json --output-dir out/work/primary"

    -- GLM primary data
    "out/work/primary/glm/singlefactor/no-outliers/Nucleotides/Muscle/Ref/PvT.csv" %> \out -> do
      need ["src/R/glm/primary-pop-compare.R", "out/work/primary/opls/no-outliers/Nucleotides/Muscle/Ref/PvT.csv"]
      cmd_ "Rscript ./src/R/glm/primary-pop-compare.R"

    -- significance table primary
    "out/supp/no-outliers/primary-pop-compare-significance.xlsx" %> \out -> do
      need ["src/python/tables/primary-pop-compare.py", "out/work/primary/glm/singlefactor/no-outliers/Nucleotides/Muscle/Ref/PvT.csv"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/tables/primary-pop-compare.py --outlier no-outliers --level 0.05 --output out/supp/no-outliers/primary-pop-compare-significance.xlsx"
    "out/supp/outliers/primary-pop-compare-significance.xlsx" %> \out -> do
      need ["src/python/tables/primary-pop-compare.py", "out/work/primary/glm/singlefactor/no-outliers/Nucleotides/Muscle/Ref/PvT.csv"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/tables/primary-pop-compare.py --outlier outliers --level 0.05 --output out/supp/outliers/primary-pop-compare-significance.xlsx"

    "out/fig/no-outliers/ratios-combined.pdf" %> \out -> do
      need ["src/python/ratio-plot-combined.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/ratio-plot-combined.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --sample-sheet ./data/primary/sample-sheet.csv --lipids-normalized ./data/lipids/normalized --lipidmaps-fa ./data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --output ./out/fig/no-outliers/ratios-combined.pdf"

    -- OPLS lipids
    "out/work/lipids/opls/outliers/category/Sphingolipids/Muscle/Ref/PvT.csv" %> \out -> do
      need ["src/python/opls/lipids-pop-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/lipids-pop-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --output-dir out/work/lipids"

    -- GLM lipids
    "out/work/lipids/glm/singlefactor/no-outliers/Muscle/Ref/PvT.csv" %> \out -> do
      need ["src/R/glm/lipids-pop-compare.R", "out/work/lipids/opls/outliers/category/Sphingolipids/Muscle/Ref/PvT.csv"]
      cmd_ "Rscript ./src/R/glm/lipids-pop-compare.R"

    -- OPLS lipid cats
    "out/work/lipidcats/opls/Liver/Ref/Classes/PvT.csv" %> \out -> do
      need ["src/python/opls/lipidcats-pop-compare.py","data/lipids/normalized/positive/liver.csv","data/lipidmaps/lipidmaps-20200724.json","data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/lipidcats-pop-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --lipidmaps-fa ./data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json --output-dir out/work/lipidcats"

    -- GLM lipid cats
--     "out/work/lipids/glm/singlefactor/no-outliers/Muscle/Ref/PvT.csv" %> \out -> do
--       need ["src/R/glm/lipids-pop-compare.R", "out/work/lipids/opls/outliers/category/Sphingolipids/Muscle/Ref/PvT.csv"]
--       cmd_ "Rscript ./src/R/glm/lipids-pop-compare.R"

