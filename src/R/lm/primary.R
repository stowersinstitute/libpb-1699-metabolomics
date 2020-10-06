library(conflicted)
library(readr)
library(tidyr)
library(dplyr)
library(performance)
library(car)
library(lme4)

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

model <- lmer(Raw_mTIC ~ (1|Category) + (1|Population) + (1|Tissue) + (1|Condition), data=primary)
# anova <- Anova(model, type = "III")
print(r2_nakagawa(model))
print(r2(model))
# print(model_performance(model))
