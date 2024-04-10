# Make an animation for a single molecule specified by command line argument,
# showing density differences across fields
library(tidyverse)
library(cowplot)
library(gganimate)
library(glue)

density_difference_inpath <- commandArgs(trailingOnly = TRUE)[1]

# Formula should be in the file name
basename <- basename(density_difference_inpath)
# Remove extension
formula <- str_extract(basename, '^[^.]+')

# Path to write the animation to
outpath <- glue('density_animation/{formula}.gif')

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

# Make the animation directly with gganimate instead of saving the frames
# separately
density_difference_framelabeled <- density_difference |>
    # The no field simulation is still in there, but we don't know where to put
    # it, so just remove it
    filter(!is.na(field_value)) |>
    # Add a "frame" column with the field number (the rank of the field value)
    # Can't just use the rank function because each field value is present
    # multiple times
    mutate(frame = as.integer(factor(field_value)))
animation <- density_difference_framelabeled |>
    ggplot(aes(x = x, y = z,
               z = density_difference, fill = density_difference)) +
    geom_raster() +
    geom_contour(color = 'white', breaks = seq(-.1, .1, by = 0.01)) +
    coord_equal() +
    scale_fill_viridis_c(rescaler = function(x, ...)
            (pmin(pmax(x, -.1), .1)*10 + 1)/2
    ) +
    theme_nothing() +
    # Following these instructions to remove margins:
    # https://stackoverflow.com/questions/31254533/when-using-ggplot-in-r-how-do-i-remove-margins-surrounding-the-plot-area
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x = NULL, y = NULL) +
    # Need to do manual transitions because the tweeting screws up the contours
    transition_manual(frame)

# Render gif
renderer <- gifski_renderer(file = outpath, width = 500, height = 500, loop = TRUE)
# Get this error if I don't specify number of frames:
# gifski only supports png files
nframe = max(density_difference_framelabeled$frame)
# If I don't assign this to a variable, it opens the file upon completion
garbage_assignment <- 
    animate(animation, renderer = renderer, rewind = TRUE,
        duration = 10, nframe = nframe)
