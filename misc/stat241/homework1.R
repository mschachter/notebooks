
steady_state = function(n, p)
{
  q = 1 - p  
  T = array(0, c(5, 5))
  T[1, 2] = p; T[2, 3] = p; T[3, 4] = p; T[4, 5] = p; T[5, 1] = p;
  
  T[1, 5] = q; T[2, 1] = q; T[3, 2] = q; T[4, 3] = q; T[5, 4] = q  ;
  Tn = T
  for (k in 1:n) { Tn = Tn%*%T }
  print(T)
  print(Tn)
}
