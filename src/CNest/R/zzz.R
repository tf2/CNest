.First.lib <- function(lib,pkg) {
   library.dynam("CNest",pkg,lib)
   cat("CNest 1.0 loaded\n")
}

