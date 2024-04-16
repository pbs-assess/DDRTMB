## code to prepare `pcod2020` dataset goes here



#read in the input files that were used in the 2020 pcod assessment. 

pcod2020dat<-read.data.file("docs/old_model/pcod.dat")
usethis::use_data(pcod2020dat, overwrite = TRUE)


pcod2020ctl<-read.control.file("docs/old_model/pcod.ctl",
                              num.gears =6,
                              num.age.gears = 1,)
usethis::use_data(pcod2020ctl, overwrite = TRUE)

pcod2020pfc<-read.projection.file("docs/old_model/pcod.pfc")
usethis::use_data(pcod2020pfc, overwrite = TRUE)




