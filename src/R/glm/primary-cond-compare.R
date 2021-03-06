library("arm")
library(MASS)

pops = c('Pachon', 'Tinaja', 'Surface')
tissues = c('Brain', 'Muscle', 'Liver')
outliers = c('with-outliers')
comparisons = c('30vR','4vR','30v4')
conditions = c('30d Starved','4d Starved','Refed')
categories = c("Aminoacids","Carbohydrates_-CCM","Fattyacids","Misc._-_sec.metabolites","Nucleotides")
tests = c("CvS")

input_type = "opls"
# input_type = "zscore"

for (tissue in tissues) {
  for (outlier in outliers) {
    for (category in categories) {
      for (pop in pops) {
        for (comp in comparisons) {
          data <- read.csv(sprintf("out/work/primary/%s/%s/%s/%s/%s/%s.csv",input_type,outlier,category,tissue,pop,comp))
          data_names <- read.csv(sprintf("out/work/primary/%s/%s/%s/%s/%s/%s.csv",input_type,outlier,category,tissue,pop,comp),check.names=FALSE)
          kegg = data[1,]
          data <- data[-1,]
          rownames(data) <- NULL
          for (k in 2:length(names(data))) {
            data[,k] <- as.double(as.character(data[,k]))
          }

          names(data)[1] <- "Condition"

          # single factor model
          coefs <- c()
          for (cpd in colnames(data)[-1]) {
      #       https://www.r-bloggers.com/changing-the-variable-inside-an-r-formula/
            if (comp != '30v4'){
              singlefac <- bayesglm(as.formula(sprintf("relevel(Condition,\"Refed\") ~ %s+0", cpd)), data = data, family = binomial())
            } else {
              singlefac <- bayesglm(as.formula(sprintf("relevel(Condition,\"4d Starved\") ~ %s+0", cpd)), data = data, family = binomial())
            }
            singlefac.summary <- summary(singlefac)
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
          dir.create(sprintf("out/work/primary/glm/singlefactor/%s/%s/%s/%s",outlier,category,tissue,pop), showWarnings = FALSE, recursive = TRUE, mode = "0755")
          write.csv(singlefac.data,sprintf("out/work/primary/glm/singlefactor/%s/%s/%s/%s/%s.csv",outlier,category,tissue,pop,comp))
        }
      }
    }
  }
}
