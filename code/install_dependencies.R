## Install the following packages

dependencies <- c("glmnet", "monomvn", "grpreg", "remotes", "bayesm", "mcclust")
for (d in dependencies) {
  if (!require(d)) install.packages(d)
}

remotes::install_github("cran/EBglmnet")
remotes::install_github("cran/GreedyEPL")
install.packages("https://cran.r-project.org/src/contrib/Archive/effectFusion/effectFusion_1.1.3.tar.gz",
                 repos = NULL, type = "source")