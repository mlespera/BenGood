## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE---------------------------------------------------------
#  library(BenGood)
#  
#  # Create dataset
#  set.seed(123)  #optional, set random number generator seed for replicable results
#  x <- rnorm(100)
#  
#  # Obtain first digit frequencies
#  freq <- firstdigitsfreq(x, 1)
#  
#  CVMStats(freq) # Cramer-von Mises test statistics
#  LR.dec(freq) # Alternative: decreasing probabilities
#  LR.genben(freq) # Alternative: Generalized Benford
#  LR.rod(freq) # Alternative: Rodriguez's generalized Benford
#  LR.mat.mult(freq) # Alternative: general multinomial, i.e. non-negative probabilities
#  
#  # 95% Simultaneous confidence intervals: for this example, Goodman and Sison
#  myCI <- SimultConf(freq, c('Goodman', 'Sis'), 0.05)
#  
#  # Plot Goodman and Sison simultaneous CIs
#  plotCI(myCI$Lower, myCI$Upper)
#  

## ------------------------------------------------------------------------
library(BenGood)

# Create dataset
set.seed(123) #optional, set random number generator seed for replicable results
x <- rnorm(100)

# Obtain first digit frequencies
freq <- firstdigitsfreq(x, 1)
# or equivalently, freqdigits(extractdigits(data))

## ------------------------------------------------------------------------
CVMStats(freq)  # Cramer-von Mises test statistics

## ------------------------------------------------------------------------
LR.dec(freq) # Alternative: decreasing probabilities
LR.genben(freq) # Alternative: Generalized Benford
LR.rod(freq) # Alternative: Rodriguez's generalized Benford
LR.mat.mult(freq) # Alternative: general multinomial, i.e. non-negative probabilities

## ------------------------------------------------------------------------
# 95% Simultaneous confidence intervals: for this example, Goodman and Sison
# Lower and upper bound matrices
myCI <- SimultConf(freq, c('Goodman', 'Sis'), 0.05)

## ----fig.align = 'center', fig.dim = c(5,5)------------------------------
# Plot Goodman and Sison simultaneous CIs
# Note: set argument legloc = TRUE if your machine supports locator 
# to choose the location of your legend
plotCI(myCI$Lower, myCI$Upper)

