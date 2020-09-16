library("arm")
library(MASS)

pops = c('Pachon', 'Tinaja', 'Surface')
tissues = c('Brain', 'Muscle', 'Liver')
outliers = c('mit-Ausreißern','kein-Ausreißern')
comparisons = c('PvS','TvS','PvT')
conditions = c('30d','4d','Ref')
categories = c("Aminoacids","Carbohydrates_-CCM","Fattyacids","Misc._-_sec.metabolites","Nucleotides")

input_type = "opls"
# input_type = "zscore"

for (tissue in tissues) {
  for (outlier in outliers) {
    for (category in categories) {
      for (cond in conditions) {
        for (comp in comparisons) {
#           https://stackoverflow.com/questions/3411201/specifying-column-names-in-a-data-frame-changes-spaces-to
          data <- read.csv(sprintf("out/work/primary/%s/%s/%s/%s/%s/%s.csv",input_type,outlier,category,tissue,cond,comp))
          data_names <- read.csv(sprintf("out/work/primary/%s/%s/%s/%s/%s/%s.csv",input_type,outlier,category,tissue,cond,comp),check.names=FALSE)
          kegg = data[1,]
          data <- data[-1,]
          rownames(data) <- NULL
          for (k in 2:length(names(data))) {
            data[,k] <- as.double(as.character(data[,k]))
          }

          names(data)[1] <- "Population"


          # single factor model
          coefs <- c()
          for (cpd in colnames(data)[-1]) {
      #       https://www.r-bloggers.com/changing-the-variable-inside-an-r-formula/
            if (comp != 'PvT'){
              singlefac <- bayesglm(as.formula(sprintf("relevel(Population,\"Surface\") ~ %s+0", cpd)), data = data, family = binomial())
            } else {
              singlefac <- bayesglm(as.formula(sprintf("relevel(Population,\"Tinaja\") ~ %s+0", cpd)), data = data, family = binomial())
            }
            singlefac.summary <- summary(singlefac)
            coefs <- c(coefs,list(singlefac.summary$coefficients))
          }
          singlefac.data <- do.call(rbind.data.frame, coefs)
          q <- p.adjust(singlefac.data[,"Pr(>|z|)"],"BH")
          singlefac.data$q <- q
          kegg <- as.vector(kegg)
          kegg <- kegg[-1]
          singlefac.data$KEGG <- t(kegg)
#           print(colnames(data_names)[-1])
          rownames(singlefac.data) <- colnames(data_names)[-1]
          dir.create(sprintf("out/work/primary/glm/singlefactor/%s/%s/%s/%s",outlier,category,tissue,cond), showWarnings = FALSE, recursive = TRUE, mode = "0755")
          write.csv(singlefac.data,sprintf("out/work/primary/glm/singlefactor/%s/%s/%s/%s/%s.csv",outlier,category,tissue,cond,comp))
#           print(sprintf("wrote out/work/primary/glm/singlefactor/%s/%s/%s/%s/%s.csv",outlier,category,tissue,cond,comp))
      }
      }
    }
  }
}
