library("arm")
library(MASS)

pops = c('Pachon', 'Tinaja', 'Surface')
tissues = c('Brain', 'Muscle', 'Liver')
outliers = c('outlier','no-outlier')
polarities = c('positive','negative')
# conditions = c('4d Starved', '30d Starved', 'Refed')
conditions = c('30d','4d','Ref')
cattypes = c('Categories','Classes')
comparisons = c('PvS','TvS','PvT')

input_type = "opls"
outlier <- "with-outliers"

# population-based comparison
for (tissue in tissues) {
  for (condition in conditions) {
    for (comparison in comparisons) {
      data <- read.csv(sprintf("out/work/lipidcats/%s/%s/%s/%s/fas/%s.csv",input_type,outlier,tissue,condition,comparison))
      names(data)[1] <- "Population"
      data$Population = factor(data$Population)

      # single factor model
      coefs <- c()
      for (cpd in colnames(data)[-1]) {
    #       https://www.r-bloggers.com/changing-the-variable-inside-an-r-formula/
        if (comparison != "PvT") {
          singlefac <- bayesglm(as.formula(sprintf("relevel(Population,\"Surface\") ~ %s+0", cpd)), data = data, family = binomial())
        } else {
          singlefac <- bayesglm(as.formula(sprintf("relevel(Population,\"Tinaja\") ~ %s+0", cpd)), data = data, family = binomial())
        }
        singlefac.summary <- summary(singlefac)
    #       https://stackoverflow.com/questions/26508519/how-to-add-elements-to-a-list-in-r-loop
    #       https://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time-o1
        coefs <- c(coefs,list(singlefac.summary$coefficients))
      }
    #     https://stackoverflow.com/questions/4227223/convert-a-list-to-a-data-frame
      singlefac.data <- do.call(rbind.data.frame, coefs)
    #     https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
      q <- p.adjust(singlefac.data[,"Pr(>|z|)"],"BH")
      singlefac.data$q <- q
#         print(singlefac.data)
      dir.create(sprintf("out/work/lipidcats/glm/%s/%s/%s/fas",outlier,tissue,condition), showWarnings = FALSE, recursive = TRUE, mode = "0755")#//
      write.csv(singlefac.data,sprintf("out/work/lipidcats/glm/%s/%s/%s/fas/%s.csv",outlier,tissue,condition,comparison))
    }
  }
}
