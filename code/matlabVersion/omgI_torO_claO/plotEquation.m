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

omgI = 314.15926;%input angular velocity [1/s] 

torO = 62.5;    %brief outpout torque [Nm]

claO = 20000;    %brief output clamping [N] 

mustar=mu/sin(delta0);

phi =1.872316189229377;

rI =0.030687090662822;

rO =0.076717726657055;

x=0.1:0.1:20.1;
y=0.1:1:150.1;
[v0,omgO]=meshgrid(x,y);

%% Equation 2,214 2,217
LOin=computeLOin(E,A,rO,omgO,v0);
K=computeK(E,A,v0,mstar);
PhiO=computePhiO(torO,rO,LOin,K,mustar);

Equation_1=rI*omgI - v0 .* (1 + ( (LOin- K) .* exp(mustar*PhiO) + K ) / (E*A));
Equation_2=claO - (E*A-mstar*v0.^2) .* (2*LOin*phi + (LOin-K) .* ((exp(mustar*PhiO)-1)/mustar - PhiO)) / (2*tan(delta0)*E*A)  +  mstar*v0.^2*phi/tan(delta0);

%% Plot


F0=zeros(size(v0));
surf(v0,omgO,Equation_1);
xlabel('v0'),ylabel('torO'),
view([-31,62]),hold on,
axis([0 20 0 100 -0.001 0.001])
surf(v0,torO,Equation_2/1e9),
caxis([-0.001,0.001])
title('Equation\_2 devided by 1e9'),
surf(v0,torO,F0),
shading interp,


%% Plot zero contour
figure,

con=[-0.001,0,0.001];

contour(v0,torO,Equation_1,con),
hold on,
contour(v0,torO,Equation_2,con),
legend('Equation_1','Equation_2')
title('Contours of Equations at height 0'),

[v0_start,torO_start]=ginput(1);

v0_start,
torO_start,

