function [omgO,v0,PhiI,PhiO,LOin]=computeUnknowns(E,A,mstar,delta0,mu,phi,rI,rO,omgI,torI,claO)




if (torI~=0)
   mustar=mu/sin(delta0);
      
 
   torO= - torI*rO/rI;   
   %omgO_start=rI*omgI/rO;
   v0_start=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
   omgO_start=rI*omgI/rO-v0_start*torO/(rO*rO*E*A);
   
   %mgO_start=1.256211094881533e+02;%1.256637040000000e+02;
   
   x0=[omgO_start,v0_start];
   
   options=optimset('TolX',1e-12,'MaxFunEvals',1000,'MaxIter',1000,'PlotFcns',@optimplotx); 
   x=fsolve(@equation_movement,x0,options);
   
   
   omgO=x(1);
   v0=x(2);
   
   LOin=computeLOin(E,A,rO,omgO,v0);
   PhiO=computePhiO(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
   PhiI=PhiO;
    
else
    
   omgO=rI*omgI/rO;
   v0=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
   PhiI=0;
   PhiO=0;
   LOin=computeLOin(E,A,rO,omgO,v0);
   
end
    
