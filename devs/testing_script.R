
# Data files have already been loaded into the package using this script

#document and build packages (these are also buttons in Rstudio)
devtools::document()
devtools::load_all()

# read in pcod data

source("devs/load-models.R")

pcod2020dat<-read.data.file("data-raw/pcod.dat")
pcod2020ctl<-read.control.file("data-raw/pcod.ctl",
                              num.gears =6,
                              num.age.gears = 1,)
pcod2020pfc<-read.projection.file("data-raw/pcod.pfc")

# TODO: read in rep file and mcmc files - put in package space

# XXX

#save these to a data folder
usethis::use_data_raw()
