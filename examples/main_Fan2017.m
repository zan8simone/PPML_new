clearvars;
close all; 
clc;
addpath('..\2d_basic_general');
addpath('..\shared_functions');


%%%%%%%%%%%%
%
% Metagrating deflector. 
% From "Large-Angle, Multifunctional Metagratings Based on Freeform
% Multimode Geometries", 
% David Sell, Jianji Yang, Sage Doshay, Rui Yang, and Jonathan A. Fan
% Nano Letters 2017 17 (6), 3752-3757
% DOI: 10.1021/acs.nanolett.7b01082 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
%
%  Simone Zanotto, 2026

%% creating the lattice

maskimp = imread('Fan_SI.png');
mask = maskimp(:,1:404,1);
mask = mask';

a1 = [1087,0];
a2 = [0,525];

extra.plotflag = 1;
extra.plotfignum = 40;
extra.plotmn = 1;
trunc.scheme = 'parall';
trunc.rad1 = 10;     
trunc.rad2 = 4;     

lattice  = create_lattice(a1,a2,size(mask,1),size(mask,2),trunc,extra); clearvars extra;
%% creating the pattern

thr = 0.5;
epsmat = 3.5^2;
epsmask.xx = (ones(size(mask)) + imbinarize(mask,thr)*(epsmat-1));
epsmask.yy = epsmask.xx;
epsmask.z  = epsmask.xx;
epsmask.xy = zeros(size(mask));
epsmask.yx = zeros(size(mask));

plotpars.plotflag = 1;
plotpars.plotfignum = 50;
plotpars.N1t = 1; plotpars.N2t = 2;
plot_epsmask(lattice,epsmask,plotpars);           clearvars plotpars;

ftpatt = create_fteps(lattice,epsmask);

%% creating the layers and solving RCWA

theta  = 0.01;
lambda = 1050;
phi    = 0;     % angle in degrees 


struc.fts    = {          ftpatt                };
struc.d      = [3000       300            3000  ]; % nm
struc.sigma  = [       0              0         ];    

struc.epssup = 1;                                
struc.epssub = 1.46^2;


globpars.k0   = 2*pi/lambda;         % wavevector in nm ^-1
inlayer = 'sub';
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
incpars.Esp = [1;0]; % [1;0] s-pol (TE)    [0;1] p-pol (TM)   
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

%% Calculating the diffraction Jones matrices
extra.plotflag = 1;
extra.plotfignum = 21;
difdata = DJM_2d_basic_general(wm,lattice,struc,globpars,extra);

%% Calculating the fields

incpars.Iinc    = 1;    
incpars.inlayer = inlayer; 
incpars.Esp     = [0,1];

fieldpars.type = 'a1_a2';
fieldpars.N1 = 100;
fieldpars.N2 = 50;
fieldpars.l = 2;
fieldpars.zz = 150;

fielddata = field_2d_basic_general(wm,lattice,struc,globpars,incpars,fieldpars);

figure(30)
subplot(121)
h = pcolor(fielddata.X,fielddata.Y,abs(fielddata.Ex).^2 + abs(fielddata.Ey).^2 + abs(fielddata.Ez).^2); hold on;
h.EdgeColor = 'none';
axis equal;
ax = gca; ax.CLim = [0 7]; ax.YDir = 'reverse';

edgeIm = edge(mask');
h = imagesc([0 a1(1)],[0 a2(2)],cat(3, ones(size(edgeIm)), zeros(size(edgeIm)), zeros(size(edgeIm))));
set(h, 'AlphaData', edgeIm);


fieldpars.type = 'a1_z';
fieldpars.N1 = 100;
fieldpars.N2 = 300;
fieldpars.f = 0.5;

fielddata = field_2d_basic_general(wm,lattice,struc,globpars,incpars,fieldpars);

subplot(122)
h = pcolor(fielddata.X,fielddata.Z,abs(fielddata.Ex).^2 + abs(fielddata.Ey).^2 + abs(fielddata.Ez).^2); hold on;
h.EdgeColor = 'none';
%quiver(fielddata.X,fielddata.Z,real(fielddata.Ex) ,real(fielddata.Ez) )
plot([0 a1(1)],[1 1]*struc.d(1),'r');
plot([0 a1(1)],[1 1]*(struc.d(1) + struc.d(2)),'r');
axis equal; ax = gca; ax.YDir = 'reverse'; ax.CLim = [0 7];