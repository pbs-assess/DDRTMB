

#document and build packages (these are also buttons in Rstudio)

devtools::document()
devtools::load_all()





# read in pcod data

source("devs/load-models.R")



pcod2020dat<-read.data.file("notes/old_model/pcod.dat")
pcod2020ctl<-read.control.file("notes/old_model/pcod.ctl",
                              num.gears =6,
                              num.age.gears = 1,)
pcod2020pfc<-read.projection.file("docs/old_model/pcod.pfc")


#save these to a data folder
usethis::use_data_raw()