function [omgO,v0,PhiI,PhiO,LOin]=computeUnknowns(E,A,mstar,delta0,mu,phi,rI,rO,omgI,torO,claO)




if (torO~=0)
   mustar=mu/sin(delta0);
 
   
   %omgO_start=rI*omgI/rO;
   v0_start=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
   omgO_start=rI*omgI/rO-v0_start*torO/(rO*rO*E*A);
   
   x0=[omgO_start,v0_start];
   
%   options=optimset('Display','iter');
 
   x=fsolve(@equation_movement,x0);
   
   
   omgO=x(1);
   v0=x(2);
   %K=computeK(E,A,v0,mstar);
   LOin=computeLOin(E,A,rO,omgO,v0);
   PhiO=computePhiO(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
   PhiI=PhiO;
   
   
%    K-computeK(E,A,v0,mstar)
%    rI*omgI-v0*(1+((LOin-K)*exp(mustar*PhiO)+K)/(E*A))
%    torO-rO*(LOin-K)*(exp(mustar*PhiO)-1)
%    claO-(E*A-mstar*v0^2)*(2*LOin*phi+(LOin-K)((exp(mustar*PhiO)-1)/mustar-PhiO))/(2*tan(delta0)*E*A)+mstar*v0^2*phi/tan(delta0)
    
else
    
   omgO=rI*omgI/rO;
   v0=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
   PhiI=0;
   PhiO=0;
   LOin=computeLOin(E,A,rO,omgO,v0);
   
end
    
