function [omgI,v0,PhiI,PhiO,LIin]=computeUnknowns(E,A,mstar,delta0,mu,phi,rI,rO,torI,omgO,claI)




if (torI~=0)
    
   mustar=mu/sin(delta0);
   
   torO=-rO*torI/rI;
   %omgI_start=rO*omgO/rI;
   v0_start=( sqrt( (-claI*tan(delta0)+E*A*(pi-phi))^2 + 4*E*A*mstar*(rO*omgO*(pi-phi))^2 ) - E*A*(pi-phi) + claI*tan(delta0) ) / (2*mstar*rO*omgO*(pi-phi));
   omgI_start=(v0_start*torO/(E*A*rO)+rO*omgO )/rI;
%    omgO_start=125;
%    v0_start=9.0;
   
   x0=[omgI_start,v0_start];
   
%   options=optimset('Display','iter');
 
   x=fsolve(@equation_movement,x0);
   
   
   omgI=x(1);
   v0=x(2);
   %K=computeK(E,A,v0,mstar);
   LIin=computeLIin(E,A,rI,omgI,v0);
   PhiI=computePhiI(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
   PhiO=PhiI;
   
   
%    K-computeK(E,A,v0,mstar)
%    rI*omgI-v0*(1+((LOin-K)*exp(mustar*PhiO)+K)/(E*A))
%    torO-rO*(LOin-K)*(exp(mustar*PhiO)-1)
%    claO-(E*A-mstar*v0^2)*(2*LOin*phi+(LOin-K)((exp(mustar*PhiO)-1)/mustar-PhiO))/(2*tan(delta0)*E*A)+mstar*v0^2*phi/tan(delta0)
    
else
    
   omgI=rO*omgO/r;
   v0=( sqrt( (E*A)^2 + 4*E*A*mstar*(rO*omgO)^2 ) - E*A ) / (2*mstar*rO*omgO);
   PhiI=0;
   PhiO=0;
   %K=computeK(E,A,v0,mstar);
   LIin=computeLIin(E,A,rI,omgI,v0);
   
end
    
