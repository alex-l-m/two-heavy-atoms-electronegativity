library(tidyverse)
library(robustbase)
library(cowplot)
library(ggdark)
library(tidymodels)
theme_set(dark_mode(theme_cowplot(font_size = 16)))
library(glue)
library(ggrepel)

iqa_energy <- read_csv('group_iqa.csv.gz', col_types = cols(
  combination_id = col_character(),
  formula = col_character(),
  donor_or_acceptor = col_character(),
  symbol = col_character(),
  total_iqa_energy = col_double()
)) |>
    rename(energy = total_iqa_energy) |>
    select(combination_id, donor_or_acceptor, symbol, energy)

bader_charge <- read_csv('bader_charge_group.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    total_bader_charge = col_double()
)) |>
    rename(charge = total_bader_charge) |>
    select(combination_id, donor_or_acceptor, symbol, charge)

charge_energy <- left_join(bader_charge, iqa_energy,
                           by = c('combination_id', 'donor_or_acceptor', 'symbol'))

# field_energies_plot <- charge_energy |>
#     group_by(symbol) |>
#     mutate(zerod_energy = energy - min(energy, na.rm = TRUE)) |>
#     ungroup() |>
#     ggplot(aes(x = charge, y = zerod_energy, color = mol_id)) +
#     facet_wrap(vars(symbol)) +
#     geom_point() +
#     xlab('Bader charge') +
#     ylab('IQA energy - min (eV)') +
#     theme(axis.text.x = element_text(angle = 45)) +
#     geom_smooth(method = lmrob, formula = y ~ x + I(x^2), se = FALSE)
# ggsave('field_energies.png', field_energies_plot, width = unit(11.5, 'in'), height = unit(4.76, 'in'))


point_energy_derivatives <- read_csv('iqa_energy_derivatives.csv.gz', col_types = cols(
    combination_id = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    charge = col_double(),
    derivative = col_double()
))

simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer()
))

neutral_nofield_combinations <- simulation_table |>
    # Filter for the combinations with no field (and no net charge)
    # For the Q-Chem simulations, "field_value" really means "CDFT constrained
    # charge". So no field applied means no CDFT constrained charge, which is
    # represented as a missing value
    filter(is.na(field_value) & molecule_charge == 0) |>
    select(combination_id)
    
# Make a table of the element of the donor and acceptor in each formula
# There is no reason this has to be in this specific script, and it could be
# useful elsewhere, so I might move it to a separate script
donor_acceptor_elements <- neutral_nofield_combinations |>
    inner_join(bader_charge, by = 'combination_id') |>
    left_join(select(simulation_table, combination_id, formula), by = c('combination_id')) |>
    select(formula, donor_or_acceptor, symbol) |>
    # Pivot so there's one row per formula
    pivot_wider(names_from = donor_or_acceptor, values_from = symbol,
                names_prefix = 'element_')

nofield_derivatives <- simulation_table |>
    inner_join(neutral_nofield_combinations, by = 'combination_id') |>
    inner_join(point_energy_derivatives, by = 'combination_id') |>
    select(formula, donor_or_acceptor, symbol, charge, derivative) |>
    mutate(`dEnergy/dCharge` = derivative,
           source = glue('{formula}:{symbol}:{donor_or_acceptor}'))

nofield_energies <- simulation_table |>
    inner_join(neutral_nofield_combinations, by = 'combination_id') |>
    inner_join(charge_energy, by = 'combination_id') |>
    select(formula, donor_or_acceptor, symbol, charge, energy) |>
    mutate(source = glue('{formula}:{symbol}:{donor_or_acceptor}')) |>
    # Only include the formulas in the derivatives table. This is a hack but I
    # don't remember how I filtered it
    inner_join(select(nofield_derivatives, formula), by = 'formula') |>
    # Another hack: everything is duplicated, no doubt due to a bug
    # upstream, but for now just deduplicating
    distinct()

write_csv(nofield_energies, 'iqa_nofield_energies.csv')



smoothed_iqa_energy <- read_csv('iqa_smoothed_energy.csv.gz', col_types = cols(
    formula = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    charge = col_double(),
    energy = col_double(),
    derivative = col_double()
)) |>
    mutate(`dEnergy/dCharge` = derivative,
           source = glue('{formula}:{symbol}:{donor_or_acceptor}'))

# Plot of the IQA atomic energies
energy_plot <- smoothed_iqa_energy |>
    # Charges between -1 and 1
    filter(-1 < charge & charge < 1) |>
    ggplot(aes(x = charge, y = energy, group = source, label = formula)) +
    facet_wrap(vars(symbol), scales = 'free_y', ncol = 4) +
    geom_line() +
    geom_point(data = nofield_energies) +
    geom_text_repel(data = nofield_energies, max.overlaps = 3, color = 'white') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('iqa_energy.png', energy_plot, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Plot of the derivatives
# Probably these element labels should just be in the original tables?
smoothed_iqa_energy_elementlabeled <- smoothed_iqa_energy |>
    left_join(donor_acceptor_elements, by = 'formula') |>
    mutate(other_group_symbol = ifelse(donor_or_acceptor == 'donor',
                                       element_acceptor,
                                       element_donor))
nofield_derivatives_elementlabeled <- nofield_derivatives |>
    left_join(donor_acceptor_elements, by = 'formula') |>
    mutate(other_group_symbol = ifelse(donor_or_acceptor == 'donor',
                                       element_acceptor,
                                       element_donor))

derivatives_plot <- smoothed_iqa_energy_elementlabeled |>
    filter(-1 < charge & charge < 1) |>
    ggplot(aes(x = charge, y = `dEnergy/dCharge`,
               group = source, color = other_group_symbol)) +
    facet_wrap(vars(symbol), ncol = 4) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90)) +
    # Vertical line to indicate the location of zero charge
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_point(data = nofield_derivatives_elementlabeled)
ggsave('iqa_energy_derivatives.png', derivatives_plot, width = unit(11.5, 'in'), height = unit(4.76, 'in'))
