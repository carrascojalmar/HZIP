## ## for creating a package
.package.Name <- "HZIP"

##.First.lib <- function(lib,pkg){
##  library.dynam(.package.Name,
##                pkg,
##                lib)
##}

##.Last.lib <- function(libpath){
##  library.dynam.unload(chname="HZIP",libpath=libpath)
##}

.onAttach <- function(...){
  packageStartupMessage("Classes and Methods for R originally developed in the")
  packageStartupMessage("Complex Statistical Modeling Laboratory (CoSMo)")
  packageStartupMessage("Department of Statistics")
  packageStartupMessage("Federal University of Bahia, Brazil (2025),")
  packageStartupMessage("by and under the direction of Jalmar M. F. Carrasco,")
  packageStartupMessage("with contributions from collaborators and students.")
  packageStartupMessage("Main functions: hzip, rHZIP.")
}

.onUnload <- function(libpath){
  library.dynam.unload("HZIP",libpath=libpath)
}
