clear all;close all;format long; 

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

omgO = 1.256466683335845e+02;%ouput angular velocity [1/s] 

claI = - 1.783169997826813e+04;    %brief output clamping [N] 

torI = - 25;    %brief outpout torque [Nm]

mustar=mu/sin(delta0);

phi =1.872316189229377;

rI =0.030687090662822;

rO =0.076717726657055;

torO=-torI*rO/rI;

x=9:0.1:10.1;
y=310:1:317;
[v0,omgI]=meshgrid(x,y);

%% Equation 2,214 2,217
% LOin=computeLOin(E,A,rO,omgO,v0);
% K=computeK(E,A,v0,mstar);
% PhiO=computePhiO(torO,rO,LOin,K,mustar);
PhiI=computePhiI(rI,omgI,rO,omgO,E,A,mstar,v0,mustar);

Equation_1=rI*omgI-v0*torO/(rO*E*A)-rO*omgO;
Equation_2=2*claI*tan(delta0) - (E*A-mstar*v0.^2) .* (2*(pi-phi)*(rI*omgI./v0-1) + torI./(-mustar*rI*E*A) - (rI*omgI./v0 - E*A./(E*A-mstar*v0.^2)).*PhiI)  +  2*mstar*v0.^2*(pi-phi);

%% Plot


F0=zeros(size(v0));
surf(v0,omgI,Equation_1);
xlabel('v0'),ylabel('torO'),
view([-31,62]),hold on,
axis([0 20 0 100 -0.001 0.001])
surf(v0,omgI,Equation_2/1e9),
caxis([-0.001,0.001])
title('Equation\_2 devided by 1e9'),
surf(v0,ogmI,F0),
shading interp,


%% Plot zero contour
figure,

con=[-0.001,0,0.001];

contour(v0,omgI,Equation_1,con),
hold on,
contour(v0,omgI,Equation_2,con),
legend('Equation_1','Equation_2')
title('Contours of Equations at height 0'),

[v0_start,omgI_start]=ginput(1),


