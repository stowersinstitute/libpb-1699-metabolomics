-- https://stackoverflow.com/questions/17287080/multi-input-multi-output-compilers-with-shake
-- https://stackoverflow.com/questions/53566429/shake-how-to-set-an-environment-variable-upon-invocation
-- https://mmhaskell.com/blog/2017/3/6/making-sense-of-multiple-monads
-- https://github.com/ndmitchell/shake/issues/170
-- https://github.com/ndmitchell/shake/issues/430#issuecomment-196054985
-- https://stackoverflow.com/questions/35938956/haskell-shake-special-rule-for-building-directories


import Development.Shake hiding ((*>))
import Development.Shake.FilePath

main = shakeArgs shakeOptions $ do
    want [-- Web API
           "out/generated/cts/client.py"
          -- ID mapping
         , "out/work/ids/kegg-to-chebi.json"
          -- Figs
         , "out/fig/with-outliers/pca-global-pop.pdf"
         , "out/fig/with-outliers/pca-categorized-primary.pdf"
         , "out/fig/with-outliers/pca-categorized-lipids.pdf"
         , "out/fig/with-outliers/pca-categorized-legend.pdf"
--          , "out/fig/with-outliers/orotic-acid-plot.pdf"
         , "out/fig/with-outliers/sugar-phosphate-plot.pdf"
         , "out/fig/with-outliers/sugar-phosphate-heatmap.pdf"
         , "out/fig/with-outliers/metabolites-of-interest.pdf"
         , "out/fig/with-outliers/primary-shared-starvation-response-30vR.pdf"
         , "out/fig/with-outliers/venn-diagrams-cross-pop.pdf"
         , "out/table/primary/starvation-response/with-outliers/CvS/30vR.tex"
         , "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Ref/PvT.csv" -- OPLS primary cross-pop
         , "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/Ref/PvT.csv" -- GLM primary cross-pop
         , "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Pachon/30vR.csv" -- OPLS primary starvation response
         , "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/CvS/30vR.csv" -- GLM primary cross-pop
         , "out/work/primary/merged-mtic.csv" -- Merged primary mTIC values
         , "out/work/lipids/merged-lipids.csv" -- Merged primary mTIC values
         , "out/work/lipids/merged-cross-pop.csv" -- Merged lipids significance
         , "out/supp/with-outliers/primary-pop-compare-significance.xlsx"
         , "out/fig/with-outliers/ratios-combined.pdf"
         , "out/work/lipids/opls/with-outliers/category/Sphingolipids/Muscle/Ref/PvT.csv" -- Lipid OPLS
         , "out/work/lipids/glm/singlefactor/with-outliers/Muscle/Ref/PvT.csv" -- Lipid GLM
         , "out/work/lipids/glm/singlefactor/with-outliers/Muscle/Pachon/30vR.csv" -- Lipid GLM (compare conditions)
         , "out/work/lipids/opls/with-outliers/category/Sphingolipids/Muscle/Pachon/30vR.csv" -- Lipid OPLS for condition comparison
         , "out/work/lipidcats/opls/with-outliers/Liver/Ref/Classes/PvT.csv" -- Lipid cats/classes OPLS
         , "out/work/lipidcats/glm/with-outliers/Muscle/Ref/Categories/PvT.csv" -- Lipid cats/classes GLM
         ]

    let primary_src = "data/primary/metabolomics-corrected.csv"
    let lipid_src = "data/lipids/normalized/positive/liver.csv"

    -- Web API for CTS
    "out/generated/cts/client.py" %> \out -> do
      need ["src/yaml/cts.yaml"]
      cmd_ "pipenv run swagger_codegen generate src/yaml/cts.yaml out/generated/cts"

    -- Get mapping to ChEBI ids
    "out/work/ids/kegg-to-chebi.json" %> \out -> do
      need [primary_src]
      cmd_ (AddEnv "PYTHONPATH" "./out/generated") "pipenv run python3 ./src/python/kegg-to-chebi.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --out out/work/ids/kegg-to-chebi.json"

    -- global pca
    "out/fig/with-outliers/pca-global-pop.pdf" %> \out -> do
      need [primary_src, "src/python/pca-global.py"]
--       PIPENV_VENV_IN_PROJECT
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-global.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet ./data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --exclude-outlier False --out-dir ./out/fig/with-outliers"

    -- categorized pca for primary
    "out/fig/with-outliers/pca-categorized-primary.pdf" %> \out -> do
      need [primary_src, "src/python/pca-categorized-primary.py"]
--       PIPENV_VENV_IN_PROJECT
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-primary.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --sample-sheet data/primary/sample-sheet.csv --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --output ./out/fig/with-outliers/pca-categorized-primary.pdf"

    -- categorized pca plot for lipids
    "out/fig/with-outliers/pca-categorized-lipids.pdf" %> \out -> do
      need [lipid_src, "src/python/pca-categorized-lipids.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-lipids.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --exclude-outlier False --output ./out/fig/with-outliers/pca-categorized-lipids.pdf"

    -- legend for categorized pca plots
    "out/fig/with-outliers/pca-categorized-legend.pdf" %> \out -> do
      need [primary_src, "src/python/pca-categorized-legend.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/pca-categorized-legend.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --sample-sheet data/primary/sample-sheet.csv --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --output ./out/fig/with-outliers/pca-categorized-legend.pdf"

    -- orotic acid plot
--     "out/fig/with-outliers/orotic-acid-plot.pdf" %> \out -> do
--       need ["src/python/orotic-acid-plot.py"]
--       cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/orotic-acid-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --output ./out/fig/with-outliers/orotic-acid-plot.pdf"

    -- sugar phosphate plot
    "out/fig/with-outliers/sugar-phosphate-plot.pdf" %> \out -> do
      need [primary_src, "src/python/sugar-phosphate-plot.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/sugar-phosphate-plot.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --output ./out/fig/with-outliers/sugar-phosphate-plot.pdf"

    -- sugar phosphate heatmap
    "out/fig/with-outliers/sugar-phosphate-heatmap.pdf" %> \out -> do
      need [primary_src, "src/python/sugar-phosphate-heatmap.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/sugar-phosphate-heatmap.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --level 0.05 --output ./out/fig/with-outliers/sugar-phosphate-heatmap.pdf --output-cbar ./out/fig/with-outliers/sugar-phosphate-heatmap-cbar.pdf"

    -- metabolites of interest
    "out/fig/with-outliers/metabolites-of-interest.pdf" %> \out -> do
      need [primary_src, "src/python/metabolites-of-interest.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/metabolites-of-interest.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --output ./out/fig/with-outliers/metabolites-of-interest.pdf --output-dir-extra ./out/fig/with-outliers/"

    -- OPLS primary cross-pop
    "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Ref/PvT.csv" %> \out -> do
      need [primary_src, "src/python/opls/primary-pop-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/primary-pop-compare.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet data/primary/sample-sheet.csv --hmdb ./data/hmdb/hmdb.json --compounds ./data/kegg/compounds.json --output-dir out/work/primary"

    -- GLM primary cross-pop
    "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/Ref/PvT.csv" %> \out -> do
      need [primary_src, "src/R/glm/primary-pop-compare.R", "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Ref/PvT.csv"]
      cmd_ "Rscript ./src/R/glm/primary-pop-compare.R"

    -- OPLS primary starvation response
    "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Pachon/30vR.csv" %> \out -> do
      need [primary_src, "src/python/opls/primary-starvation-response.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/primary-starvation-response.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet data/primary/sample-sheet.csv --hmdb ./data/hmdb/hmdb.json --compounds ./data/kegg/compounds.json --output-dir out/work/primary"

    -- GLM primary condition compare
    "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/Pachon/30vR.csv" %> \out -> do
      need ["src/R/glm/primary-cond-compare.R", "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Ref/PvT.csv"]
      cmd_ "Rscript ./src/R/glm/primary-cond-compare.R"

    -- GLM primary starvation response
    "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/CvS/30vR.csv" %> \out -> do
      need ["src/R/glm/primary-starvation-response.R", "out/work/primary/opls/with-outliers/Nucleotides/Muscle/Pachon/30vR.csv"]
      cmd_ "Rscript ./src/R/glm/primary-starvation-response.R"

    -- Table of tidy data for primary mTIC values
    "out/work/primary/merged-mtic.csv" %> \out -> do
      need [primary_src, "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/CvS/30vR.csv", "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/Pachon/30vR.csv", "src/python/primary-merged.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/primary-merged.py --astyanax ./data/primary/metabolomics-corrected.csv --sample-sheet ./data/primary/sample-sheet.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --exclude-outlier False --kegg-to-chebi out/work/ids/kegg-to-chebi.json --out-mtic ./out/work/primary/merged-mtic.csv --out-cross-pop ./out/work/primary/merged-cross-pop.csv --out-starvation-resp ./out/work/primary/merged-starvation-resp.csv"

    -- Table of tidy data for lipid mTIC values
    "out/work/lipids/merged-lipids.csv" %> \out -> do
      need [lipid_src, "src/python/lipids-merged.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/lipids-merged.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --lipidmaps-fa ./data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json --out ./out/work/lipids/merged-lipids.csv"

    -- Table of tidy data for lipid significance values
    "out/work/lipids/merged-cross-pop.csv" %> \out -> do
      need ["out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/CvS/30vR.csv", "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/Pachon/30vR.csv", "src/python/lipids-stats.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/lipids-stats.py --exclude-outlier False --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --lipidmaps-fa ./data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json --out-cross-pop ./out/work/lipids/merged-cross-pop.csv --out-cross-cond ./out/work/lipids/merged-cross-cond.csv"

    -- conserved metabolites in starvation resistance
    "out/fig/with-outliers/primary-shared-starvation-response-30vR.pdf" %> \out -> do
      need ["src/python/primary-shared-starvation-response.py", "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/CvS/30vR.csv"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/primary-shared-starvation-response.py --exclude-outlier False --output-dir ./out/fig"

    -- table for conserved metabolites in starvation resistance
    "out/table/primary/starvation-response/with-outliers/CvS/30vR.tex" %> \out -> do
      need ["src/python/primary-shared-starvation-response.py", "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/CvS/30vR.csv"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/primary-shared-starvation-response-tables.py --exclude-outlier False --out-dir ./out/table/primary/starvation-response/with-outliers/CvS"

    -- conserved metabolites in starvation resistance
    "out/fig/with-outliers/venn-diagrams-cross-pop.pdf" %> \out -> do
      need [primary_src, "src/python/venn/primary-pop-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/venn/primary-pop-compare.py --exclude-outlier False --level 0.05 --output-dir ./out/fig"

    -- significance table primary
    "out/supp/with-outliers/primary-pop-compare-significance.xlsx" %> \out -> do
      need ["src/python/tables/primary-pop-compare.py", "out/work/primary/glm/singlefactor/with-outliers/Nucleotides/Muscle/Ref/PvT.csv"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/tables/primary-pop-compare.py --outlier with-outliers --level 0.05 --output out/supp/with-outliers/primary-pop-compare-significance.xlsx"

    "out/fig/with-outliers/ratios-combined.pdf" %> \out -> do
      need [primary_src, "src/python/ratio-plot-combined.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/ratio-plot-combined.py --astyanax ./data/primary/metabolomics-corrected.csv --compounds ./data/kegg/compounds.json --hmdb ./data/hmdb/hmdb.json --sample-sheet ./data/primary/sample-sheet.csv --lipids-normalized ./data/lipids/normalized --lipidmaps-fa ./data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --output ./out/fig/with-outliers/ratios-combined.pdf --output-legend ./out/fig/with-outliers/ratios-combined-legend.pdf"

    -- OPLS lipids compare pop
    "out/work/lipids/opls/with-outliers/category/Sphingolipids/Muscle/Ref/PvT.csv" %> \out -> do
      need [primary_src, "src/python/opls/lipids-pop-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/lipids-pop-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --output-dir out/work/lipids"

    -- GLM lipids compare pop
    "out/work/lipids/glm/singlefactor/with-outliers/Muscle/Ref/PvT.csv" %> \out -> do
      need ["src/R/glm/lipids-pop-compare.R", "out/work/lipids/opls/with-outliers/category/Sphingolipids/Muscle/Ref/PvT.csv"]
      cmd_ "Rscript ./src/R/glm/lipids-pop-compare.R"

    -- OPLS lipids compare condition
    "out/work/lipids/opls/with-outliers/category/Sphingolipids/Muscle/Pachon/30vR.csv" %> \out -> do
      need [primary_src, "src/python/opls/lipids-cond-compare.py"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/lipids-cond-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --output-dir out/work/lipids"

    -- GLM lipids compare condition
    "out/work/lipids/glm/singlefactor/with-outliers/Muscle/Pachon/30vR.csv" %> \out -> do
      need ["src/R/glm/lipids-cond-compare.R", "out/work/lipids/opls/with-outliers/category/Sphingolipids/Muscle/Pachon/30vR.csv"]
      cmd_ "Rscript ./src/R/glm/lipids-cond-compare.R"

    -- OPLS lipid cats
    "out/work/lipidcats/opls/with-outliers/Liver/Ref/Classes/PvT.csv" %> \out -> do
      need ["src/python/opls/lipidcats-pop-compare.py","data/lipids/normalized/positive/liver.csv","data/lipidmaps/lipidmaps-20200724.json","data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json"]
      cmd_ (AddEnv "PYTHONPATH" "./src/python") "pipenv run python3 ./src/python/opls/lipidcats-pop-compare.py --lipids-normalized ./data/lipids/normalized --lipidmaps-json ./data/lipidmaps/lipidmaps-20200724.json --lipidmaps-fa ./data/lipidmaps/lipidmaps-20200724-sat-unsat-fas.json --output-dir out/work/lipidcats"

    -- GLM lipid cats
    "out/work/lipidcats/glm/with-outliers/Muscle/Ref/Categories/PvT.csv" %> \out -> do
      need ["src/R/glm/lipid-cats-pop-compare.R", "out/work/lipidcats/opls/with-outliers/Liver/Ref/Classes/PvT.csv"]
      cmd_ "Rscript ./src/R/glm/lipid-cats-pop-compare.R"

