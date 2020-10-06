library(conflicted)
library(readr)
library(tidyr)
library(dplyr)
library(car)

conflict_prefer("filter", "dplyr")

primary <- read_csv("out/work/primary/merged-mtic.csv") %>%
  arrange(KEGG) %>%
  filter(Outlier == FALSE) %>%
  select(!Outlier)
print(primary)
