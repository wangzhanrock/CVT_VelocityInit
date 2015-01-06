function [rI,rO,phi,v0,PhiI,PhiO,LOin]=computeUnknowns(E,A,mstar,delta0,mu,lR_tilde,dA,d_align,omgI,omgO,claI,claO)




% if (omgO/omgI~=rI/rO)
 
    mustar=mu/sin(delta0);
    
  rI_start=(0.5*lR_tilde-dA)/pi;
  rO_start=(0.5*lR_tilde-dA)/pi;
  phi_start=0.5*pi;
 % v0_start=( sqrt( (claO*tan(delta0)+E*A*phi_start)^2 + 4*E*A*mstar*(rI_start*omgI*phi_start)^2 ) - E*A*phi_start - claO*tan(delta0) ) / (2*mstar*rI_start*omgI*phi_start);
%  rO_start=0.076717726657055;    
  % phi_start=1.072758155568609e+02;    
  %v0_start=9.635874561660197;
%  rI_start=0.030687090662822;  
%  rI_start=(0.5*lR_tilde-sqrt(dA^2+d_align^2))/pi;
 % rO_start=(0.5*lR_tilde-sqrt(dA^2+d_align^2))/pi;

  v0_start=( sqrt( (-claI*tan(delta0)+E*A*(pi-phi_start))^2 + 4*E*A*mstar*(rI_start*omgI*(pi-phi_start))^2 ) - E*A*(pi-phi_start) + claI*tan(delta0) ) / (2*mstar*rI_start*omgI*(pi-phi_start));
    
    x0=[rI_start,rO_start,phi_start,v0_start];
    
    %options=optimset('PlotFcns',@optimplotx);
  %  options=optimset('PlotFcns',@optimplotfval);
    options=optimset('TolX',1e-12,'MaxFunEvals',1000,'MaxIter',1000,'PlotFcns',@optimplotx);
    
    x=fsolve(@equation_movement,x0,options);
    
    
    rI=x(1);
    rO=x(2);
    phi=x(3);
    v0=x(4);
    

    LOin=computeLOin(E,A,rO,omgO,v0);
    PhiO=computePhiO(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);
    PhiI=PhiO;
    
    
    
%else
    
%    v0=( sqrt( (claO*tan(delta0)+E*A*phi)^2 + 4*E*A*mstar*(rI*omgI*phi)^2 ) - E*A*phi - claO*tan(delta0) ) / (2*mstar*rI*omgI*phi);
%   PhiI=0;
%    PhiO=0;
%    LOin=computeLOin(E,A,rO,omgO,v0);
    
end
