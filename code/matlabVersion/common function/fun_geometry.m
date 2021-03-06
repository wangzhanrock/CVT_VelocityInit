function [phi,rI,rO]=fun_geometry(ir,dA,lR_tilde,d_align)

equation_geometry=@(x)([x(2)-ir*x(2)-cos(x(1))*dA;2*sqrt(sin(x(1))^2*dA^2+d_align^2)+2*x(2)*(pi-x(1))+2*ir*x(2)*x(1)-lR_tilde]);
%equation_geometry=@(phi,rI)([rI-ir*rI-cos(phi)*dA;2*sqrt(sin(phi)^2*dA^2+d_align^2)+2*rI*(pi-phi)+2*ir*rI*phi-lR_tilde]);

 
 
rI_start=(0.5*lR_tilde-dA)/pi;
phi_start=0.5*pi;
%options=optimset('Display','iter');
x=fsolve(equation_geometry,[phi_start; rI_start]);

phi=x(1);
rI=x(2);
rO=ir*rI;