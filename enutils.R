box::use(tikzDevice[tikz])
box::use(glue[glue])
box::use(readr[read_lines])
box::use(ggplot2[ggsave])
box::use(systemfonts[match_font])
# Path to the font I'm using so that the latex output matches
font_path <- match_font('sans')$path
# Latex is very picky, needs to have the slash at the end or I get an error
font_dir <- glue('{dirname(font_path)}/')
font_base <- basename(font_path)
bold_font_base <- basename(match_font('sans', bold = TRUE)$path)
italic_font_base <- basename(match_font('sans', italic = TRUE)$path)
bold_italic_font_base <- basename(match_font('sans', bold = TRUE, italic = TRUE)$path)
# Definitions of custom commands that I use in my theory section and may also
# want to use in plot labels
preamble_lines <- read_lines('manuscript/preamble.tex')
# Lines to get fonts how I want them
font_lines <- c(
    '\\usepackage{fontspec}',
    glue('\\setmainfont[Path={font_dir}, BoldFont={{{bold_font_base}}}, ItalicFont={{{italic_font_base}}}, BoldItalicFont={{{bold_italic_font_base}}}]{{{font_base}}}'),
    glue('\\setsansfont[Path={font_dir}, BoldFont={{{bold_font_base}}}, ItalicFont={{{italic_font_base}}}, BoldItalicFont={{{bold_italic_font_base}}}]{{{font_base}}}')
)
# unique is just so I can run this and an interactive session without errors
# the second time
options(
    tikzXelatexPackages = unique(c(getOption('tikzXelatexPackages'), '\\usepackage{amsmath}', preamble_lines, font_lines))
)
# Closure providing arguments to the tikz device
device_function <- function(...) {
    tikz(file = ..., standAlone = TRUE,
    engine = 'xetex')
}

texsave <- function(outbase, plot, width, height)
{
    texfile <- glue('{outbase}.tex')
    ggsave(texfile, plot,
           width = width, height = height,
            device = device_function)
    system2('xelatex', args = texfile)
    system2('pdftoppm', args = c('-r', '600', '-png', '-singlefile', glue('{outbase}.pdf'), outbase))
}
