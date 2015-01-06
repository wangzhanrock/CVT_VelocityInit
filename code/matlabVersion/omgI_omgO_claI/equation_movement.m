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

%% input condition
ir   = 2.5;

omgI = 314.15926;%input angular velocity [1/s] 

omgO = 1.256466683335845e+02;

claI = -1.783169997826725e+04    %brief output clamping [N] 


mustar=mu/sin(delta0);

[phi,rI,rO]=fun_geometry(ir,dA,lR_tilde,d_align);


%% Unknowns
torI=x(1);
v0=x(2);
%torO= - torI*rO/rI;
%LIin=computeLIin(E,A,rI,omgI,v0);
%K=computeK(E,A,v0,mstar);
PhiI=computePhiI(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);

%% Equation 2,214 FCI

Equation(1)=rI*omgI + v0*torI/(rI*E*A)-rO*omgO;
Equation(2)=2*claI*tan(delta0) + (E*A-mstar*v0^2) * (2*(pi-phi)*(rI*omgI/v0-1) + torI/(-mustar*rI*E*A) - (rI*omgI/v0 - E*A/(E*A-mstar*v0^2))*PhiI)  -  2*mstar*v0^2*(pi-phi);

end

