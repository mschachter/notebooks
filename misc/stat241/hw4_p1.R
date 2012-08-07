

em_iterate = function(pi0, niter, y)
{
  pi = pi0
  
  for (k in 1:niter)
  {
    #do E-step
    x2 = ((0.25*pi) / (0.5 + 0.25*pi))*y[1]
    
    #do M-step
    #pi_new = (y[2] + y[3]) / (y[2] + y[3] - x2 - y[4])
    pi_new = (x2 + y[4]) / (x2 + y[4] - y[2] - y[3])
    
    print(sprintf('x2=%0.2f, old_pi=%0.3f, new_pi=%0.3f', k, x2, pi, pi_new))
    
    pi = pi_new
  }
  
  
  
  
  
  
}

