library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
theme_set(dark_mode(theme_cowplot()))

smoothed_energy <- read_csv("smoothed_energy.csv.gz", col_types = cols(
    formula = col_character(),
    symbol_donor = col_character(),
    symbol_acceptor = col_character(),
    charge_transfer = col_double(),
    energy = col_double(),
    derivative = col_double()
))

this_theme <- 
    theme(
        # x axis text is too crowded, rotate it
        axis.text.x = element_text(angle = 90)
    )

energy_comparison_comparison <- smoothed_energy |>
    ggplot(mapping = aes(x = charge_transfer, y = energy)) +
    facet_wrap(~ formula, scales = "free", nrow = 2) +
    geom_smooth(method = lmrob, formula = y ~ x + I(x^2), se = FALSE) +
    geom_line() +
    this_theme
ggsave("energy_comparison_comparison.png", energy_comparison_comparison, width = unit(11.5, "in"), height = unit(4.76, "in"))

energy_derivatives_comparison <- smoothed_energy |>
    ggplot(aes(x = charge_transfer, y = derivative)) +
    facet_wrap(vars(formula), scales = "free", nrow = 2) +
    geom_line() +
    this_theme
ggsave("energy_derivatives_comparison.png", energy_derivatives_comparison, width = unit(11.5, "in"), height = unit(4.76, "in"))
