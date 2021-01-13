.First.lib <- function(lib,pkg) {
   library.dynam("Rbin",pkg,lib)
   cat("roots 0.1-1 loaded\n")
}

