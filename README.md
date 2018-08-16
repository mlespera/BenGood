# BenGood
R package for testing goodness-of-fit of Benford's Law for first digits.  See Lesperance M, Reed WJ, Stephens MA, Tsao C, Wilton B
(2016) Assessing conformance with Benford's Law: goodness-of-fit tests and simultaneous confidence intervals. PLoS one, 11(3).

To install from GitHub, first

    install.packages("devtools")

    devtools::install_github("mlespera/BenGood", build_vignettes = TRUE)
    
If the vignette does not install, try re-installing (possibly several times),

    devtools::install_github("mlespera/BenGood", build_vignettes = TRUE, force = TRUE)
