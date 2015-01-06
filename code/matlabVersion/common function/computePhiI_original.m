function PhiI=computePhiI(torI,rI,LIin,K,mustar)
 
  PhiI = - log( 1 - torI/(rI*(LIin-K)))/mustar;
 
 
 end