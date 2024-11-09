# Use gganimate to make animations of slices of cube files
library(tidyverse)
library(glue)
library(gganimate)
library(cowplot)
library(ggdark)

slices <- read_csv('slices.csv.gz', col_types = cols(
    cube_id = col_character(),
    i = col_double(),
    j = col_double(),
    k = col_double(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    density = col_double(),
    simulation_id = col_character(),
    field_number = col_integer(),
    structure_id = col_character(),
    cube_file_type = col_character(),
    cube_file_path = col_character()
))

# Create the output directory if it doesn't already exist
outdir <- 'slice_animations'
if (!dir.exists(outdir))
{
    dir.create(outdir)
}

# Cube file types to animate
target_types <- c('cube_file_path',
    'hartree_pot_path',
    'density_derivatives_path',
    'dhartree_path',
    'pot_file_path')

# Loop over structures
slices_list <- slices |>
    # Filter for just the kinds of cube file that I want to animate
    filter(cube_file_type %in% target_types) |>
    group_by(structure_id, cube_file_type) |>
    group_split()
for (this_slices in slices_list)
{
    animation <- this_slices |>
        # Filter out anything too far above the 95% quantile
        # This helps when plotting densities, which have spikes that throw off
        # the scale
        group_by(simulation_id) |>
        filter(density < 3 * quantile(density, 0.95)) |>
        ungroup() |>
        ggplot(aes(x, y, color = density)) +
        geom_point() +
        scale_color_viridis_c() +
        # Manual frames with no tweening
        transition_manual(field_number) +
        coord_fixed() +
        dark_mode(theme_nothing()) +
        # Include a legend on the bottom, and put the title back in the theme
        theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5)) +
        # Actually don't include the title because it makes it take forever for
        # some reason
        # Don't display name of variable in the legend, because I'm using the
        # name "density" for various stuff that isn't really the density
        guides(color = guide_colorbar(title = NULL))

     
    structure_id <- this_slices$structure_id[1]
    cube_type <- this_slices$cube_file_type[1]
    anim_save(glue('{outdir}/slices_{structure_id}_{cube_type}.gif'),
              animation, duration = 1, width = 500)
}
