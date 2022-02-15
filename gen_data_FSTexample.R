library(CompQuadForm)
library(SKAT)
library(Matrix)
library(MASS)
library(mvtnorm)
# Depends R (>= 2.10), CompQuadForm, SKAT, Matrix, MASS, mvtnorm
library(FSTpackage)
# Data example for FSTest
# The dataset contains outcome variable Y, covariate X, genotype data G, 
# functional scores Z and gene-set ID for each variable GeneSetID.
data(FST.example)

gen_data_new3 = function(label.size, m1, m0, s1, s0){
  # extract Y, X, G from the sample dataset
  x = FST.example$X
  y = FST.example$Y
  g = FST.example$G
  # nn = 5088, whcih is the number of record in example
  nn = length(y)
  # Lable ids
  id.t = sample(1:nn, label.size)
  # I think the functional scores Z is not the score we need
  # so I keeped this part from origional function
  # Generate surrogates
  s = rep(NA, nn)
  s[y == 1] <- rnorm(n = sum(y), mean = m1, sd = s1)
  s[y == 0] <- rnorm(n = sum(1-y), mean = m0, sd = s0)
  # Return data
  my_dat <- cbind(y, x, g, s)
  colnames(my_dat) <- c("Y", "X", "G", "S")  # this line need revise
  # I think X and G from example data set is ready for dataframe, but error occur
  return(list(my_dat = my_dat, id.t = id.t))
}

# test --------------------------------------------------------------------

n_lab = 100 # consider labeled data 100, 200, 500 for now
m1 = 2
m0 = 1
s1 = 0.5
s0 = 0.6
set.seed(001)

dat.obj = gen_data_new3(label.size = n_lab, m1, m0, s1, s0)
dat = dat.obj$my_dat
id.t = dat.obj$id.t
