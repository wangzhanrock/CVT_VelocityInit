function Equation=equation_movement(x)

%% paramerter
Ne = 432;
mE = 1.858e-3;
rhor = 7.7e3;
Nr = 9;
hrl = 1.85e-4;
rr = 1.044e-1;
wr = 9.75e-3;
hr = Nr * hrl;
lR_tilde = 2 * pi * (rr + hr / 2);



E = 1.85e11;  %brief belt Young's modulus [N/m^2]

A = hr * 2. * wr;  % brief belt cross-sectional area [m^2]

mstar = (Ne * mE + rhor * lR_tilde * hr * 2. * wr) / lR_tilde;  %brief belt mass distribution [kg/m]

mu = 0.11;

d_align = -0.0001060006182;

dA = 1.55e-1;

delta0 = 11 * pi / 180; %brief pulley basic half wedge angle [-] 

mustar=mu/sin(delta0);

%% input condition

omgI = 314.15926;%input angular velocity [1/s] 

omgO = 1.256466683335845e+02;

claI= -1.783169997826813e+04;

claO = 20000;    %brief output clamping [N] 


%%unknows

rI=x(1);
rO=x(2);
phi=x(3);
v0=x(4);


Equation(1)=rI-rO-cos(phi)*dA;    % @rI,rO,phi   #:dA
Equation(2)=2*sqrt(sin(phi)^2*dA^2+d_align^2)+2*rI*(pi-phi)+2*rO*phi-lR_tilde;  %@rI,rO,phi, #dA,d_align,lR_tilde

Equation(3)=2*claO*tan(delta0) - (E*A-mstar*v0^2)*(2*phi*(rO*omgO/v0 -1) + (rI*omgI-rO*omgO)/(mustar*v0)-(rO*omgO/v0 - E*A/(E*A-mstar*v0^2) )*log((rI*omgI*(E*A-mstar*v0^2)-E*A*v0)/(rO*omgO*(E*A-mstar*v0^2)-E*A*v0))/mustar) + 2*mstar*v0^2*phi;    %@v0,phi,rI,rO  #claO,delta0,E,A,mstar,mustar,

Equation(4)=2*claI*tan(delta0) + (E*A-mstar*v0^2)*(2*(pi-phi)*(rI*omgI/v0-1) + (rO*omgO-rI*omgI)/(-mustar*v0) - (rI*omgI/v0 - E*A/(E*A-mstar*v0^2))*log((rO*omgO*(E*A-mstar*v0^2)-E*A*v0)/(rI*omgI*(E*A-mstar*v0^2)-E*A*v0))/(-mustar))  -  2*mstar*v0^2*(pi-phi);
