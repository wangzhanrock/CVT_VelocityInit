function PhiI=computePhiI(rI,omgI,rO,omgO,E,A,mstar,v0,mustar)
 
  %PhiI=log( (rO*omgO*(E*A-mstar*v0^2)-E*A*v0)/(rI*omgI*(E*A-mstar*v0^2)-E*A*v0) )/(-mustar);
 PhiI=log( (rO*omgO*(E*A-mstar*v0.^2)-E*A*v0)./(rI*omgI.*(E*A-mstar*v0.^2)-E*A*v0) )/(-mustar);
 
 end
