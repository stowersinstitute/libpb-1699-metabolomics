library("arm")
library(MASS)

pops = c('Pachon', 'Tinaja', 'Surface')
tissues = c('Brain', 'Muscle', 'Liver')
outliers = c('with-outliers')
comparisons = c('PvS','TvS','PvT')
conditions = c('30d','4d','Ref')

input_type = "opls"
# input_type = "zscore"

for (tissue in tissues) {
  for (outlier in outliers) {
    for (cond in conditions) {
      for (comp in comparisons) {
#           https://stackoverflow.com/questions/3411201/specifying-column-names-in-a-data-frame-changes-spaces-to
        data <- read.csv(sprintf("out/work/lipids/%s/%s/%s/%s/%s.csv",input_type,outlier,tissue,cond,comp))
        data_names <- read.csv(sprintf("out/work/lipids/%s/%s/%s/%s/%s.csv",input_type,outlier,tissue,cond,comp),check.names=FALSE)
        rownames(data) <- NULL
        for (k in 2:length(names(data))) {
          data[,k] <- as.double(as.character(data[,k]))
        }

        names(data)[1] <- "Population"

        # all-factor model
#           m <- model.matrix(~ .+0, data[-1])
#           if (comp != '30v4'){
#             allfac <- bayesglm(relevel(Population,"Surface") ~ m+0, data = data, family = binomial())
#           } else {
#             allfac <- bayesglm(relevel(Population,"Tinaja") ~ m+0, data = data, family = binomial())
#           }
#           glm.summary <- summary(allfac)
#
#           dir.create(sprintf("/tmp/primary/glm/allfactor/%s/%s/%s/%s",outlier,category,tissue,comp), showWarnings = FALSE, recursive = TRUE, mode = "0755")
#           write.csv(glm.summary$coefficients,sprintf("/tmp/primary/glm/allfactor/%s/%s/%s/%s/model.csv",outlier,category,tissue,comp))

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
        rownames(singlefac.data) <- colnames(data_names)[-1]
        dir.create(sprintf("out/work/lipids/glm/singlefactor/%s/%s/%s",outlier,tissue,cond), showWarnings = FALSE, recursive = TRUE, mode = "0755")
        write.csv(singlefac.data,sprintf("out/work/lipids/glm/singlefactor/%s/%s/%s/%s.csv",outlier,tissue,cond,comp))
#         print(sprintf("wrote /tmp/lipids/glm/singlefactor/%s/%s/%s/%s.csv",outlier,tissue,cond,comp))
      }
    }
  }
}
