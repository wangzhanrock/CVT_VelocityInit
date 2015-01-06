function [torI,v0,PhiI,PhiO,LOin]=computeUnknowns(E,A,mstar,delta0,mu,phi,rI,rO,omgI,omgO,claI)




if (omgO/omgI~=rI/rO)
    
    mustar=mu/sin(delta0);
    
    v0_start=( sqrt( (-claI*tan(delta0)+E*A*(pi-phi))^2 + 4*E*A*mstar*(rI*omgI*(pi-phi))^2 ) - E*A*(pi-phi) + claI*tan(delta0) ) / (2*mstar*rI*omgI*(pi-phi));
    torI_start=-rI*E*A*(rI*omgI-rO*omgO)/v0_start;
    
    x0=[torI_start,v0_start];
    
    %  options=optimset('Display','iter');
    
    x=fsolve(@equation_movement,x0);
    
    
    torI=x(1);
    v0=x(2);
    %torO=-torI*rO/rI;
    %K=computeK(E,A,v0,mstar);
    LOin=computeLOin(E,A,rO,omgO,v0);
    PhiO=computePhiO(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
    PhiI=PhiO;
    
    
else
    
    
    v0=( sqrt( (claI*tan(delta0)+E*A*(pi-phi))^2 + 4*E*A*mstar*(rO*omgO*(pi-phi))^2 ) - E*A*(pi-phi) - claI*tan(delta0) ) / (2*mstar*rO*omgO*(pi-phi));
    PhiI=0;
    PhiO=0;
    LOin=computeLOin(E,A,rO,omgO,v0);
    
end
