function [omgO,v0,PhiI,PhiO,LIin]=computeUnknowns(E,A,mstar,delta0,mu,phi,rI,rO,omgI,torI,claI)




if (torI~=0)
    
   mustar=mu/sin(delta0);
   
   torO=-rO*torI/rI;
   %omgI_start=rO*omgO/rI;

   v0_start=( sqrt( (-claI*tan(delta0)+E*A*(pi-phi))^2 + 4*E*A*mstar*(rI*omgI*(pi-phi))^2 ) - E*A*(pi-phi) + claI*tan(delta0) ) / (2*mstar*rI*omgI*(pi-phi));
   omgO_start=rI*omgI/rO-v0_start*torO/(rI*rO*E*A);
   
   x0=[omgO_start,v0_start];
   
%   options=optimset('Display','iter');
 
   x=fsolve(@equation_movement,x0);
   
   
   omgO=x(1);
   v0=x(2);
   
   LIin=computeLIin(E,A,rI,omgI,v0);
   PhiI=computePhiI(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
   PhiO=PhiI;
   
    
else
    
   omgO=rI*omgI/r;
   v0=( sqrt( (claI*tan(delta0)+E*A*(pi-phi))^2 + 4*E*A*mstar*(rI*omgI*(pi-phi))^2 ) - E*A*(pi-phi) - claI*tan(delta0) ) / (2*mstar*rI*omgI*(pi-phi));
   PhiI=0;
   PhiO=0;
   LIin=computeLIin(E,A,rI,omgI,v0);
   
end
    
