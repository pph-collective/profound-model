# run this file to install all the dependencies for profound
to_install <- c(
    "dplyr",
    "tictoc",
    "openxlsx",
    "rstudioapi")

for (i in to_install) {
  message(paste("looking for ", i))
  if (!requireNamespace(i)) {
    message(paste("     installing", i))
    install.packages(i)
  }
}
