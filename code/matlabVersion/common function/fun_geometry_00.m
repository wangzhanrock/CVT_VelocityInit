function [phi,rI,rO]=fun_geometry_00(ir,dA,lR_tilde,d_align)

f=@(x)([x(2)-ir*x(2)-cos(x(1))*dA;2*sqrt(sin(x(1))^2*dA^2+d_align^2)+2*x(2)*(pi-x(1))+2*ir*x(2)*x(1)-lR_tilde]);
%f=@(phi,rI)([rI-ir*rI-cos(phi)*dA;2*sqrt(sin(phi)^2*dA^2+d_align^2)+2*rI*(pi-phi)+2*ir*rI*phi-lR_tilde]);

rI_start=(0.5*lR_tilde-dA)/pi;
phi_start=0.5*pi;
options=optimset('Display','iter')
x=fsolve(f,[phi_start; rI_start],options);

phi=x(1);
rI=x(2);
rO=ir*rI;


