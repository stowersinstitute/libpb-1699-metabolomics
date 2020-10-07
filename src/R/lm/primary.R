library(conflicted)
library(readr)
library(tidyr)
library(dplyr)
library(performance)
library(lme4)
library(car)

conflict_prefer("filter", "dplyr")

primary <- read_csv("out/work/primary/merged-mtic.csv",
                    col_types=cols(
                      Category = col_factor(),
                      Population = col_factor(),
                      Tissue = col_factor(),
                      Condition = col_factor()
                    )) %>%
  arrange(KEGG) %>%
  filter(Outlier == FALSE) %>%
  select(!Outlier)
print(primary)

model <- lmer(Raw_mTIC ~ Name + (1|Population) + (1|Tissue), data=primary)
print(model)
# anova <- Anova(model, type = "III")
# print(anova)
print(r2(model, by_group=TRUE))
# print(pp_check(model))
# print(coef(model))
# print(confint(model))
# print(varianceDecomposition(model))
# print(model_performance(model))
