
library(gjam)


load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

PA_data<- load_object("DOM.mat.sites.species.PA.RData")
CA_data<- load_object("DOM.mat.sites.species.PA.RData")
