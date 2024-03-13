# Make an animation for a single molecule specified by command line argument,
# showing density differences across fields
library(tidyverse)
library(cowplot)

density_difference_inpath <- commandArgs(trailingOnly = TRUE)[1]

# Formula should be in the file name
basename <- basename(density_difference_inpath)
# Remove extension
formula <- str_extract(basename, '^[^.]+')

density_difference <- read_csv(density_difference_inpath, col_types = cols(
    combination_id = col_character(),
    field_value = col_double(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    coordinate_1 = col_double(),
    coordinate_2 = col_double(),
    density = col_double(),
    density_nofield = col_double(),
    density_difference = col_double()
))
# Split into a list of tables, one for each combination, to loop over to make
# the individual frames
density_difference_list <- density_difference |>
    group_by(combination_id) |>
    group_split()

for (this_density_difference_tbl in density_difference_list)
{
    this_combination_id <- this_density_difference_tbl$combination_id[1]
    dir.create(sprintf('density_images/%s', formula), showWarnings = FALSE)
    this_density_difference_plot <- this_density_difference_tbl |>
        ggplot(aes(x = x, y = z,
                   z = density_difference, fill = density_difference)) +
        geom_raster() +
        geom_contour(color = 'white', breaks = seq(-.1, .1, by = 0.01)) +
        coord_equal() +
        scale_fill_viridis_c(rescaler = function(x, ...) (pmin(pmax(x, -.1), .1)*10 + 1)/2) +
        theme_nothing() +
        # Following these instructions to remove margins:
        # https://stackoverflow.com/questions/31254533/when-using-ggplot-in-r-how-do-i-remove-margins-surrounding-the-plot-area
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        labs(x = NULL, y = NULL)
    ggsave(sprintf('density_images/%s/%s.png', formula, this_combination_id), width = unit(5, 'in'), height = unit(5, 'in'), dpi = 100)
}
