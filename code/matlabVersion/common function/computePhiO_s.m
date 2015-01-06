%% calculate the value of PhiO 
function PhiO=computePhiO_s(torO,v0,E,A,mstar,rO,omgO,mustar)

PhiO=log( (torO*v0*(E*A-mstar*v0^2)) / ( (  rO*omgO*(E*A-mstar*v0^2) - v0*E*A) *E*A*rO )  + 1 )  /  mustar;


end

