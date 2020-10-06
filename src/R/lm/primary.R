library(conflicted)
library(readr)
library(tidyr)
library(dplyr)
# library(performance)
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
# print(primary)

model <- aov(Raw_mTIC~Category + Population + Tissue + Condition, data=primary)
anova <- Anova(model, type = "III")
print(anova)
print(r2(model))
# print(model_performance(model))
