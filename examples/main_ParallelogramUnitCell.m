clearvars;
close all; 
clc;
addpath('..\2d_basic_general');
addpath('..\shared_functions');

% Metagrating with parallelogram-shaped unit cell
%  Grating layer is a dielectric with polygonal holes
%  There is a conducting layer between superstrate and dielectric layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
%
%  Simone Zanotto, 2026

epsmat = 3.37^2;  % something like simicon
sigmaG = 6e-5;    % something like monolayer graphene
%% Initializing the lattice

a1 = [2000,0]; 
a2 = [200,1000]; 

extra.plotflag = 1;
extra.plotfignum = 40;
extra.plotmn = 1;
% trunc.scheme = 'parall';
% trunc.rad1 = 3;     
% trunc.rad2 = 1;     
trunc.scheme = 'circ';
trunc.rad = 3;   

lattice  = create_lattice(a1,a2,512,256,trunc,extra); clearvars extra;
%%
%%%%%%%%% Creating patterns and their Fourier transforms

polygon = [100,100; 1000,100; 800,600; NaN,NaN; 1600,200; 2000,400; 1600,900; 1200,400]; 

epsA.xx = epsmat; 
epsA.yy = epsmat; 
epsA.z  = epsmat; 
epsA.xy = 0; 
epsA.yx = 0;

epsB.xx = 1; 
epsB.yy = 1; 
epsB.z  = 1; 
epsB.xy = 0; 
epsB.yx = 0;

epsmask = create_epsmask_from_polygon(lattice,polygon,epsA,epsB);     

plotpars.plotfignum = 50;
plotpars.N1t = 1; plotpars.N2t = 1;
plot_epsmask(lattice,epsmask,plotpars);           clearvars plotpars;


fteps = create_fteps(lattice,epsmask);

%% Declaring layer structure and solving RCWA

lambda = 1550;

% conical incidence
theta  = 60.01; 
phi = 30;

%              superstr.    |  PG   |  subst.       
struc.d       = [5000         300      5000  ]; % 
struc.fts     = {             fteps          }; % material A 
struc.sigma   = [    sigmaG               0         ];    

struc.epssup = 1^2;                                
struc.epssub = 1.4^2;  

globpars.k0   = 2*pi/lambda;         % wavevector in nm ^-1

inlayer = 'sup';
switch inlayer
    case 'sup'  
    globpars.kparx = globpars.k0 * sqrt(struc.epssup) * sind(theta)*cosd(phi);
    globpars.kpary = globpars.k0 * sqrt(struc.epssup) * sind(theta)*sind(phi);
    case 'sub'
    globpars.kparx = globpars.k0 * sqrt(struc.epssub) * sind(theta)*cosd(phi);
    globpars.kpary = globpars.k0 * sqrt(struc.epssub) * sind(theta)*sind(phi);
end

wm = solve_2d_basic_general(lattice,struc,globpars);

%% Calculating the diffracted intensities
incpars.inlayer = inlayer;
incpars.Esp = [1;-1i]; % [1;0] s-pol   [0;1] p-pol   [1;1i] R-pol   [1;-1i] L-pol 
                      %   (following   Chipman-Tiffany-Young, CRC Press 2019)
extra.plotflag = 1;
extra.plotfignum = 20;
RTAdata = RTA_2d_basic_general(wm,lattice,struc,globpars,incpars,extra);

% writing the efficiencies for all orders alongside with the (m,n) order indexing
[RTAdata.R, lattice.mn]
[RTAdata.T, lattice.mn]

% writing the efficiencies only for open orders alongside with the (m,n) order indexing
switch inlayer
    case 'sup'
    [RTAdata.R(RTAdata.IIsup), lattice.mn(RTAdata.IIsup,:)]
    [RTAdata.T(RTAdata.IIsub), lattice.mn(RTAdata.IIsub,:)]
    case 'sub'
    [RTAdata.R(RTAdata.IIsub), lattice.mn(RTAdata.IIsub,:)]
    [RTAdata.T(RTAdata.IIsup), lattice.mn(RTAdata.IIsup,:)]
end

% writing the absorption layer-by-layer
RTAdata.A