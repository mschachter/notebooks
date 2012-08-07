
pca1 = read.table('pca1.dat')
pca2 = read.table('pca2.dat')

pca = function(ds)
{
  nrows = dim(ds)[1]
  ncols = dim(ds)[2]
  
  cmat = cov(cbind(ds[, 1], ds[, 2]))
  eset = eigen(cmat)
  e1 = eset$vectors[, 1]
  
  pcadata = list(lambda=e1, psi=diag(ncols))
  
  return(pcadata)
}


outer_sum = function(ds)
{
  nrows = dim(ds)[1]
  ncols = dim(ds)[2]
  
  opsum = mat.or.vec(ncols, ncols)
  for (k in 1:nrows) {
    x = as.numeric(ds[k, ])
    opsum = opsum + (x %o% x)
  }
  return(opsum)
}

factor_analysis = function(ds, p, lambda0, psi0, niter=10)
{
  nrows = dim(ds)[1]
  ncols = dim(ds)[2]
  
  lambda = lambda0
  psi = psi0
  dsmean = mean(ds)
  os = outer_sum(ds)
  I = diag(p)
  
  for (k in 1:niter) {
    
    #do E-step
    ia = solve((lambda %*% t(lambda)) + psi)
  
    ipsi = solve(psi)
    var_Xy = solve(I + as.numeric(lambda %*% ipsi %*% lambda))    
    #var_Xy = I - ((lambda %*% solve((lambda %o% lambda) + psi)) %*% lambda)
    
    #e_Xy = as.matrix(ds[, ] - dsmean) %*% t(lambda %*% ia)
    e_Xy = as.matrix(ds[, ] - dsmean) %*% t(var_Xy %*% lambda %*% ipsi)    
    e_XX = var_Xy + e_Xy[, ] * e_Xy[, ]
    
    #do M step for lambda
    b = as.numeric(colSums(ds * e_Xy))    
    lambda = as.numeric(b * solve(sum(e_XX)))
    
    #do M step for psi    
    psi = diag(ncols) * (diag(os - (lambda %o% b)) / nrows)    
  }
  
  fadata = list(lambda=lambda, psi=psi)
  
  return(fadata)
}


compare_pca_and_fa = function(ds)
{
  nrows = dim(ds)[1]
  ncols = dim(ds)[2]
  
  #generate random initial guesses for factor analysis
  lambda0 = runif(2)
  psi0 = diag(runif(2))
  fadata = factor_analysis(ds, 1, lambda0, psi0, niter=1000)
  
  #do PCA
  pcadata = pca(ds)

  #project onto subspaces
  pca_proj = as.matrix(ds) %*% pcadata$lambda 
  fa_proj = as.matrix(ds) %*% fadata$lambda
  
  #compute mean and variance of p(x | y) for each y
  I = diag(1)
  dsmean = mean(ds)
  ipsi = solve(fadata$psi)
  var_xy = solve(I + as.numeric(fadata$lambda %*% ipsi %*% fadata$lambda))
  std_xy = sqrt(var_xy)
  e_xy = as.matrix(ds[, ] - dsmean) %*% t(var_xy %*% fadata$lambda %*% ipsi)
  
  #compute posterior p(x | y) for each x
  p_xy = mat.or.vec(nrows, 1)
  for (k in 1:nrows) {    
    p_xy[k] = dnorm(fa_proj[k], mean=e_xy[k], sd=std_xy)
  }
  
  #make histograms of subspace projections
  pca_mean = mean(pca_proj)
  pca_std = sd(pca_proj)  
  fa_mean = mean(fa_proj)
  fa_std = sd(fa_proj)  
  par(mfrow=c(3, 1))
  hist(pca_proj, main=sprintf('PCA: mean=%0.2f +/- %0.2f', pca_mean, pca_std))
  hist(fa_proj, main=sprintf('Factor: mean=%0.2f +/- %0.2f', fa_mean, fa_std))
  
  #plot x vs p(x | y)
  plot(fa_proj, p_xy, xlab='x', ylab='p(x|y)')
  
  rdata = list(pcadata=pcadata, fadata=fadata)
  return(rdata)
}
