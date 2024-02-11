.onLoad <- function(lib, pkg){
   library.dynam("CNest", pkg, lib)
}
.onUnload <- function(libpath)
    library.dynam.unload("CNest", libpath)
