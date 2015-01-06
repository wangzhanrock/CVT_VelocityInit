clear all;close all;clc;format long;

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


%% computing

[rI,rO,phi,v0,PhiI,PhiO,LOin]=computeUnknowns(E,A,mstar,delta0,mu,lR_tilde,dA,d_align,omgI,omgO,claI,claO);

%% E,A,mstar,delta0,torO,claO,omgI,rO,rI,phi,v0,omgO,PhiO,PhiI,LOin,v0,omgO,K=computeK(E,A,v0,mstar),

torO=(rI*omgI-rO*omgO)*rO*E*A/v0;
torI=-rI*torO/rO;
LIin=computeLOin(E,A,rI,omgI,v0);
K=computeK(E,A,v0,mstar);
mustar=mu/sin(delta0);
claI=computeClaI(E,A,mstar,v0,LIin,phi,K,mustar,PhiI,delta0);


%% Output

sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'),
ir=rO/rI,

rO, rO_ref=0.076717726657055,
rO-rO_ref,

rI, rI_ref=0.030687090662822,
rI-rI_ref,

phi_degree=180/pi*phi,  phi_degree_ref=1.072758155568609e+02,
phi_degree-phi_degree_ref,

v0,v0_ref=9.635874561660197,
v0-v0_ref,

PhiO_degree=180/pi*PhiO, PhiO_degree_ref=33.732112521095750,
PhiO_degree-PhiO_degree_ref,

LOin,  LOin_ref=2.151913200452924e+03,
LOin-LOin_ref,

LIin,  LIin_ref=2.966588018431328e+03,
LIin-LIin_ref,

omgO,  omgO_ref=1.256466683335845e+02,
omgO-omgO_ref,


claI,  claI_ref= -1.783169997826813e+04
claI-claI_ref,


torI,  torI_ref=-25,
torI-torI_ref,

torO,  torO_ref=62.5,
torO-torO_ref,
sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'),
