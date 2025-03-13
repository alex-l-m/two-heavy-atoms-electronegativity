# Make a plot of all of the line densities of the atomic references
library(tidyverse)
library(glue)
library(cowplot)
library(ggdark)
theme_set(dark_mode(theme_cowplot()))
library(robustbase)
library(broom)

# Function for reading the line density output from MultiWFN
# The first three columns are coordinates. I'm not sure what the fourth column
# is, but the fifth is the density. The columns are separated by an unknown
# number of spaces. There's no header row
read_line_density <- function(inpath)
{
    read_table(inpath, col_names=c('x', 'y', 'z', 'unknown', 'density'),
               col_types = cols(
                   x = col_double(),
                   y = col_double(),
                   z = col_double(),
                   unknown = col_double(),
                   density = col_double()
                 ))
}

atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_double()
))

# Line densities of all ions
line_densities <- tibble(inpath = Sys.glob('single_atoms/*_line.txt')) |>
    mutate(ion_id = str_match(basename(inpath),
                              '^([A-Z][a-z]?_[+-]?[0-9]+)_line.txt$')[,2]) |>
    group_by(ion_id) |>
    reframe(read_line_density(inpath)) |>
    mutate(element = str_extract(ion_id, '^[A-Z][a-z]?'),
           charge = as.integer(str_extract(ion_id, '[+-]?[0-9]+$')))

# Table of neutral or cation
line_density_nonnegative <- line_densities |>
    filter(charge >= 0)
# Table of neutral or anion
line_density_nonpositive <- line_densities |>
    filter(charge <= 0)
# Combining into a separated table
line_density_sep <- bind_rows(`Neutral and cation` = line_density_nonnegative,
                              `Neutral and anion` = line_density_nonpositive,
                              .id='Ion type') |>
    mutate(`Ion type` = factor(`Ion type`,
                               levels = c('Neutral and cation',
                                          'Neutral and anion')))

elements_present <- line_density_sep |>
    distinct(element)

ordered_elements <- elements_present |>
    left_join(atomic_numbers, by=c('element'='symbol')) |>
    arrange(atomic_number) |>
    pull(element)

# Plot, with each element as a facet and ions superimposed
zmax <- 5
line_density_plt <- line_density_sep |>
    # Order the elements by atomic number
    mutate(element = factor(element, levels=ordered_elements)) |>
    filter(abs(charge) <= 2 & z < zmax & 10^-5 < density & density < 10) |>
    ggplot(aes(x=z, y=density, color=as.factor(charge))) +
    facet_grid(`Ion type` ~ element) +
    geom_line() +
    scale_y_log10() +
    scale_x_continuous(limits=c(0, zmax)) +
    labs(x = 'Distance (Angstroms)', y = 'Density (unknown unit)') +
    theme(legend.position = 'bottom')

ggsave('single_atom_lines.png', line_density_plt,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Add linear fits to the plot
line_density_plt_fits <- line_density_plt +
    geom_smooth(method=lmrob, formula=y ~ x, se=FALSE, linetype='dashed')

ggsave('single_atom_lines_fits.png', line_density_plt_fits,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Save linear fits to a table using broom
line_density_fits <- line_densities |>
    filter(abs(charge) <= 2 & z < zmax & 10^-5 < density & density < 10) |>
    group_by(ion_id, element, charge) |>
    do(tidy(lmrob(log(density) ~ z, data=.))) |>
    ungroup()
write_csv(line_density_fits, 'single_atom_lines_fits.csv')

# Make a plot of the slopes
slater_exponent_plot <- line_density_fits |>
    filter(term == 'z') |>
    ggplot(aes(x=charge, y=estimate)) +
    facet_wrap(~element) +
    geom_point()
ggsave('single_atom_slopes.png', slater_exponent_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))
