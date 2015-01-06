function claI=computeClaI(E,A,mstar,v0,LIin,phi,K,mustar,PhiI,delta0)

claI= - (E*A-mstar*v0^2) * (2*LIin*(pi-phi) + (LIin-K) * ((exp(-mustar*PhiI)-1)/(-mustar) - PhiI)) / (2*tan(delta0)*E*A)  +  mstar*v0^2*(pi-phi)/tan(delta0);
