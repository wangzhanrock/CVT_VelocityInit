%% not be used.

function Equation=equation_geometry(X)
 
 ir=1;
 dA=100;
 lR_tilde=200+40*pi;
 d_align=0;
 
 phi=X(1);
 rI=X(2);

 Equation(1)=rI-ir*rI-cos(phi)*dA;
 Equation(2)=2*sqrt(sin(phi)^2*dA^2+d_align^2)+2*rI*(pi-phi)+2*ir*rI*phi-lR_tilde;

 Equation=Equation';
  
end
