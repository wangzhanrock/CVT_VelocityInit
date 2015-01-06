function PhiO=computePhiO(rI,omgI,rO,omgO,E,A,mstar,v0,mustar)
 
  PhiO=log( (rI*omgI*(E*A-mstar*v0^2)-E*A*v0)/(rO*omgO*(E*A-mstar*v0^2)-E*A*v0) )/mustar;
 
 
 end
