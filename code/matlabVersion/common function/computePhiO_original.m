function PhiO=computePhiO(torO,rO,LOin,K,mustar)
 
  PhiO=log( torO/(rO*(LOin-K))  + 1 )/mustar;
 
 
 end
