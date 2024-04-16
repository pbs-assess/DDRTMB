## code to prepare `pcod2020` dataset goes here



#read in the input files that were used in the 2020 pcod assessment. 

pcod2020dat<-read.data.file("data-raw/pcod.dat")
usethis::use_data(pcod2020dat, overwrite = TRUE)


nomedat<-names(pcod2020dat)
pcod2020dat[nomedat[2]]




pcod2020ctl<-read.control.file("data-raw/pcod.ctl",
                              num.gears =6,
                              num.age.gears = 1,)
usethis::use_data(pcod2020ctl, overwrite = TRUE)


nomectl<-names(pcod2020ctl)

t(t(names(pcod2020ctl)))

pcod2020pfc<-read.projection.file("data-raw/pcod.pfc")
usethis::use_data(pcod2020pfc, overwrite = TRUE)


nomepfc<-names(pcod2020pfc)

