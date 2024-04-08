# Extract the "Lam" parameter from lines grepped out of the Q-Chem logs

library(tidyverse)

lam_lines <- readLines('lamvals.txt')
# Example line:
# qchem_logs/C1Cl1H3_F30_0.log: Lam  -0.1690143569386403
# Parsing regex with groups for the id and the lambda value
lam_regex <- '^qchem_logs/(.*)\\.log: Lam\\s+([\\d\\.-]+)$'
# Parse to get a table
lam_values_raw <- tibble(lines = lam_lines) |>
    transmute(combination_id = str_match(lines, lam_regex)[, 2],
              lam_raw = as.numeric(str_match(lines, lam_regex)[, 3]))
# My best guess at a unit conversion
TO_EV <- 27.211386246
lam_values <- lam_values_raw |>
    mutate(lam = lam_raw * TO_EV) |>
    select(-lam_raw)

write_csv(lam_values, 'lam_values.csv.gz')
