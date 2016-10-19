# visdom-loadshape
Load shape clustering tools that extend the capabilities of the VISDOM load data analysis R package. To install, follow these steps:

1. Install the R module visdom following instructions found here: https://github.com/convergenceda/visdom and ensure that you can load the visdom package. The upshot is that you should be able to install visdom using `devtools::install_github('convergenceda/visdom')`, but you will need to ensure that you have first installed devtools and depending on your planned usage there are other ways to do it.

2. Install the visdom-loadshape package. Since it is currently a private repo, you will need to get permission to access it. 

See the code in the vignettes directory, especially the R markdown file [example_load_shape.rmd](./vignettes/example_load_shape.rmd). Note: rmd's mix commentary and code and this one will give you a sense of what is required to run load shape clustering on a given sample of meter data.

