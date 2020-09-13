library("arm")
library(MASS)

pops = c('Pachon', 'Tinaja', 'Surface')
tissues = c('Brain', 'Muscle', 'Liver')
outliers = c('outliers','no-outliers')
comparisons = c('30vR','4vR','30v4')
conditions = c('30d Starved','4d Starved','Refed')
categories = c("Aminoacids","Carbohydrates_-CCM","Fattyacids","Misc._-_sec.metabolites","Nucleotides")
tests = c("CvS")

for (tissue in tissues) {
  for (outlier in outliers) {
    for (category in categories) {
      for (test in tests) {
        for (comp in comparisons) {
          pachon <- read.csv(sprintf("out/work/primary/opls/%s/%s/%s/Pachon/%s.csv",outlier,category,tissue,comp))
          tinaja <- read.csv(sprintf("out/work/primary/opls/%s/%s/%s/Tinaja/%s.csv",outlier,category,tissue,comp))
          surface <- read.csv(sprintf("out/work/primary/opls/%s/%s/%s/Surface/%s.csv",outlier,category,tissue,comp))
          data_names <- read.csv(sprintf("out/work/primary/opls/%s/%s/%s/Surface/%s.csv",outlier,category,tissue,comp),check.names=FALSE)
          kegg = pachon[1,]
          pachon <- pachon[-1,]
          rownames(pachon) <- NULL
          tinaja <- tinaja[-1,]
          rownames(tinaja) <- NULL
          surface <- surface[-1,]
          rownames(surface) <- NULL
          for (k in 2:length(names(pachon))) {
            pachon[,k] <- as.double(as.character(pachon[,k]))
            tinaja[,k] <- as.double(as.character(tinaja[,k]))
            surface[,k] <- as.double(as.character(surface[,k]))
          }
          if (outlier == "no-outliers") {
            if (tissue == "Liver") {
              # Pachon 30d outlier
              if (comp == "30vR" || comp == "30v4") {
                tinaja <- tinaja[-3,]
                surface <-surface[-3,]
              }
              # Tinaja R outlier
              if (comp == "30vR" || comp == "4vR") {
                pachon <- pachon[-11,]
                surface <-surface[-11,]
              }
            }
            else if (tissue == "Muscle") {
              # Pachon R outlier
              if (comp == "30vR" || comp == "4vR") {
                tinaja <- tinaja[-11,]
                surface <- surface[-11,]
              }
            }
          }

          names(pachon)[1] <- "Condition"
          names(tinaja)[1] <- "Condition"
          names(surface)[1] <- "Condition"
#           print(table(pachon$Condition))
#           print(table(tinaja$Condition))
#           print(table(surface$Condition))
          for (condition in conditions){
            stopifnot(sum(pachon$Condition == condition) == sum(tinaja$Condition == condition))
            stopifnot(sum(tinaja$Condition == condition) == sum(surface$Condition == condition))
          }
          contrasted = 0.5*(pachon+tinaja) - surface
          contrasted$Condition = factor(pachon$Condition)

          # all-factor model
          m <- model.matrix(~ .+0, contrasted[-1])
          if (comp != '30v4'){
            allfac <- bayesglm(relevel(Condition,"Refed") ~ m+0, data = contrasted, family = binomial())
          } else {
            allfac <- bayesglm(relevel(Condition,"4d Starved") ~ m+0, data = contrasted, family = binomial())
          }
          glm.summary <- summary(allfac)

          dir.create(sprintf("out/work/primary/glm/allfactor/%s/%s/%s/%s",outlier,category,tissue,test,comp), showWarnings = FALSE, recursive = TRUE, mode = "0755")
          write.csv(glm.summary$coefficients,sprintf("out/work/primary/glm/allfactor/%s/%s/%s/%s.csv",outlier,category,tissue,test,comp))

          # single factor model
          coefs <- c()
          for (cpd in colnames(contrasted)[-1]) {
      #       https://www.r-bloggers.com/changing-the-variable-inside-an-r-formula/
            if (comp != '30v4'){
              singlefac <- bayesglm(as.formula(sprintf("relevel(Condition,\"Refed\") ~ %s+0", cpd)), data = contrasted, family = binomial())
            } else {
              singlefac <- bayesglm(as.formula(sprintf("relevel(Condition,\"4d Starved\") ~ %s+0", cpd)), data = contrasted, family = binomial())
            }
            singlefac.summary <- summary(singlefac)
      #       https://stackoverflow.com/questions/26508519/how-to-add-elements-to-a-list-in-r-loop
      #       https://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time-o1
            coefs <- c(coefs,list(singlefac.summary$coefficients))
          }
          singlefac.data <- do.call(rbind.data.frame, coefs)
      #     https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
          q <- p.adjust(singlefac.data[,"Pr(>|z|)"],"BH")
          singlefac.data$q <- q
          kegg <- as.vector(kegg)
          kegg <- kegg[-1]
          singlefac.data$KEGG <- t(kegg)
          rownames(singlefac.data) <- colnames(data_names)[-1]
          dir.create(sprintf("out/work/primary/glm/singlefactor/%s/%s/%s/%s",outlier,category,tissue,test), showWarnings = FALSE, recursive = TRUE, mode = "0755")
          write.csv(singlefac.data,sprintf("out/work/primary/glm/singlefactor/%s/%s/%s/%s/%s.csv",outlier,category,tissue,test,comp))
#           print(sprintf("wrote out/work/primary/glm/singlefactor/%s/%s/%s/%s/%s.csv",outlier,category,tissue,test,comp))
        }
      }
    }
  }
}
