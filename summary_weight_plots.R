
library(ggplot2)
#setwd("/mnt/beegfs/mistis/dbystrov/Gjammod")
setwd("/Users/dariabystrova/Documents/GitHub/gjamed")


###################################################################################################

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

###################################################################################################
###Load existing tables




########combine tables
#S=1000
tab_1000<-load_object("tablecompS1000_all.rds")

#S=100
tab_100<-load_object("tablecompS100.rds")

#S=200
tab_200<-load_object("tablecompS200.rds")


#S=500
tab_500<-load_object("tablecompS500.rds")


