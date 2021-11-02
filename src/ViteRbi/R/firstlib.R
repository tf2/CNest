.onLoad <- function(lib, pkg){
   library.dynam("ViteRbi", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("ViteRbi", libpath)

