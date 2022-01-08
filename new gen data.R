## Necessary libraries.
library(SPAtest)

library(MASS)
library(SKAT)
library(mltools)
if(!require(SRAT)) {
  devtools::install_github("celehs/SRAT") 
}
library(SRAT)
library(mvtnorm)
library(dplyr)
library(matrixStats)
library(sim1000G)


## Helper functions.
logit = function(xx){log(xx/(1-xx))}
g.logit = function(xx){exp(xx)/(1+exp(xx))}
autocorr.mat <- function(p = 100, rho = 0.3){
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

# Repeat vc in dm rows. 
VTM <- function(vc, dm){
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}


## Data generation function.
gen_data_new2 <- function(nn, prev=0.2, a0_cov, 
                         m1, m0, s1, s0, label.size=100, rho = 0.3){
  
  # Fix prevalence of the outcome.
  y <- rbinom(nn, 1, prev)
  
  # Generate covariates.
  tmp.mean <- VTM(a0_cov[1, ], nn) + y * VTM(a0_cov[2, ], nn)
  num_cov <- ncol(a0_cov) #+ length(b0_gene)
  x <- mvtnorm::rmvnorm(nn, sigma = autocorr.mat(num_cov, rho)) + tmp.mean
  
  # Generate genetic markers.
  # Read the example file included in sim1000G
  examples_dir = system.file("examples", package = "sim1000G")
  vcf_file = file.path(examples_dir,"region.vcf.gz")
  # Alternatively provide a vcf file here:
  # vcf_file = "~/fs/tmp/sim4/pop1/region-chr4-312-GABRB1.vcf.gz"
  # maxNumberOfVariants : total number of variants from the area that will be kept, generally should be < 1000
  # min_maf : minimum allele frequency of the markers to keep. It should be >0.
  # max_maf : maximum allele frequency of markers to keep
  # maximumNumberOfIndividuals : how many simulated individuals to reserve space for
  vcf = readVCF( vcf_file, maxNumberOfVariants = 500  , min_maf = 0.01, max_maf = 1)
  startSimulation(vcf, totalNumberOfIndividuals = nn)
  ids = generateUnrelatedIndividuals(nsim)
  g = retrieveGenotypes(ids)
  
  # Generate surrogates.
  s = rep(NA, nn)
  s[y == 1] <- rnorm(n = sum(y), mean = m1, sd = s1)
  s[y == 0] <- rnorm(n = sum(1-y), mean = m0, sd = s0)
  ## TODO: Add similar analysis for poisson and binomial.  
  ## First step: Do poisson, then transform to log(x+1) and use current
  ## assumption of normality.
  
  # Labeled ids. 
  id.t = sample(1:nn, label.size)
  
  # Return data. 
  my_dat <- cbind(y, x, g, s)
  colnames(my_dat) <- c("Y", paste0("X", 1:ncol(x)),
                        paste0("G",1:ncol(g)),
                        "S")
  
  return(list(my_dat = data.frame(my_dat), id.t = id.t))
}

# test --------------------------------------------------------------------

n_lab <- 100
## set parameters
a = c(-1, 0.4, 0.4, 0.2, 0.2, 0.3)
a0_cov = rbind(rep(a[1], length(a)-1), a[-1])

m1 = 2
m0 = 1
s1 = 0.5
s0 = 0.6
log.s1 = log(s1)
log.s0 = log(s0)

## TODO: Next, see how surrogate changes result.  Try AUC 0.5, 0.75. 
## Current setting is 0.9.

nn = 5000
#
n_lab = 100 # consider labeled data 100, 200, 500 for now
prev = 0.2
nsim = 1000
#
n_boot = 50


# function
set.seed(1234)
dat.obj = gen_data_new2(nn, prev, a0_cov, m1, m0, s1, s0, 
                        label.size = n_lab)
dat = dat.obj$my_dat
id.t = dat.obj$id.t
