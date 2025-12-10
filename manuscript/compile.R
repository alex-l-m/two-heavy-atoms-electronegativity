library(pandoc)
library(bookdown)

pandoc_activate('3.8.3')
render_book('index.Rmd')
