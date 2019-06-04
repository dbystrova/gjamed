
setwd("~/Documents/GitHub/gjamed")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


L4GJ<-load_object("Sim_smallS_gjam.Rda")