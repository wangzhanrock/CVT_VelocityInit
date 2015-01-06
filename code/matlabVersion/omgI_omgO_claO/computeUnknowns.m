function [torO,v0,PhiI,PhiO,LOin]=computeUnknowns(E,A,mstar,delta0,mu,phi,rI,rO,omgI,omgO,claO)




if (omgO/omgI~=rI/rO)
    mustar=mu/sin(delta0);
    
    
    v0_start=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
    torO_start=(rI*omgI-rO*omgO)*rO*E*A/v0_start;
    
    x0=[torO_start,v0_start];
    
    %  options=optimset('Display','iter');
    %options=optimset('TolX',1e-12,'MaxFunEvals',1000,'MaxIter',1000,'PlotFcns',@optimplotstepsize);
    
    x=fsolve(@equation_movement,x0);
    
    
    torO=x(1);
    v0=x(2);
    LOin=computeLOin(E,A,rO,omgO,v0);
    PhiO=computePhiO(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
    PhiI=PhiO;
    
    
else
    
    v0=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
    PhiI=0;
    PhiO=0;
    LOin=computeLOin(E,A,rO,omgO,v0);
    
end
