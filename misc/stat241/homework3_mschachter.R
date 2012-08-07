
ds = read.table('lms.dat', header=FALSE, col.names=c('x1', 'x2', 'y'))

solve_normal_equation = function(ds)
{
  X = as.matrix(ds[, 1:2])
  y = as.array(ds[, 3])
  
  A = t(X) %*% X
  b = y %*% X
  
  params = solve(A, t(b))
  
  return(params)
}

covariance = function(ds)
{
  X = as.matrix(ds[, 1:2])
  y = as.array(ds[, 3])
  n = nrow(X)
    
  #xmean = colMeans(X)
  #covmat = mat.or.vec(2, 2)  
  #for (k in 1:n) {   
  #  xm = X[k, ] - xmean
  #  c = outer(xm, xm)
  #  covmat = covmat + c
  #}
  #covmat = covmat / (n-1)
  covmat = cov(X)
  
  return(covmat)  
}

cost_function = function(param1, param2)
{
  nargs = length(param1)
  sumsq = mat.or.vec(nargs, 1)
  for (k in 1:nargs) {
    
    params = c(param1[k], param2[k])
    X = as.matrix(ds[, 1:2])  
    y = as.array(ds[, 3])
    n = nrow(X)
    
    ypred = X %*% params
    err = as.vector(y) - as.vector(ypred)
    sqerr = err**2
    
    sumsq[k] = sum(sqerr)
  }
  
  
  return(sumsq)
}

contour_plot = function(ds)
{
  bestparams = solve_normal_equation(ds)
  
  xgrid = seq(from=-4, to=6, by=0.2)
  ygrid = seq(from=-4, to=4, by=0.2)
  
  z = outer(xgrid, ygrid, FUN="cost_function")

  contour(xgrid, ygrid, z)  
  points(bestparams[1], bestparams[2])
}


lms_path = function(ds, step)
{
  X = as.matrix(ds[, 1:2])
  y = as.array(ds[, 3])
  n = nrow(X)
  
  path = mat.or.vec(n+1, 2)
  
  for (k in 2:n+1) {
    yn = y[k-1]
    xn = X[k-1, ]
    params = path[k-1, ]  
    
    newparams = params + step*(yn - (params %*% xn))*xn
    
    path[k, ] = newparams
  }
  return(path)
}

lms_plot = function(ds)
{
  best_params = solve_normal_equation(ds)
  covmat = covariance(ds)
  eset = eigen(covmat)  
  
  maxval = max(eset$values)
  
  step1 = 1 / maxval
  step2 = 0.5*step1
  step3 = 0.25*step1
  
  path1 = lms_path(ds, step1)
  path2 = lms_path(ds, step2)
  path3 = lms_path(ds, step3)
  
  plot(path1[, 1], path1[, 2], type="b", col=2, pch='o')
  lines(path2[, 1], path2[, 2], type="b", col=3, pch='o')
  lines(path3[, 1], path3[, 2], type="b", col=4, pch='o')  
  points(best_params[1], best_params[2], pch='+')
}

