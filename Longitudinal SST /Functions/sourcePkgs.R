# Short script to load required packages for EXG-SST analysis.

# Note: the 'dmc' folder must be in your working directory when running this script

# Load packages
rpkgs <- c("reshape2", "Hmisc", "ggplot2", "psych", "skimr",
             "cowplot", "ggpubr", "mgcv", "broom",
             "dplyr", "plyr", "tidyr", "mice", "boot",
             "SemiPar", "readxl", "gridExtra", "grid", "lattice")

if (length(setdiff(rpkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(rpkgs, rownames(installed.packages())))
}
lapply(rpkgs, library, character.only = TRUE)


# load dmc packages
dmcpkgs <- c("msm", "loo", "hypergeo", "statmod", "rtdists",
             "pracma", "snowfall", "rlecuyer", "numDeriv",
             "vioplot", "ggplot2", "gridExtra", "mvtnorm",
             "Matrix", "Brobdingnag", "stringr","LaplacesDemon",
             "coda")
lapply(dmcpkgs, library, character.only = TRUE)


# Load the EXG3 Model
# source the DMC package from the folder
source ("dmc/dmc.R")
# load the EXG model
load_model ("EXG-SS", "exgSSprobit.R")
cores = 20
