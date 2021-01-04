.onLoad <- function(lib, pkg){
   library.dynam("Rbin", pkg, lib)
}
.onUnload <- function(libpath)
    library.dynam.unload("Rbin", libpath)
