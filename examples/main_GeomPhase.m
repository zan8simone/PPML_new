clearvars;
close all; 
clc;
addpath('..\2d_basic_general');
addpath('..\shared_functions');

% Geometric phase birefringent grating

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
%
%  Simone Zanotto, 2026

%% Initializing the lattice

a1 = [6000,0]; 
a2 = [0,6000]; % NOTE: since, later on, the structure will be defined invariant with respect to y, 
               % any value of the form [0, ...] will ultimately yield the same results
               % (with a suitable remapping of diffraction orders)

extra.plotflag = 1;
extra.plotfignum = 40;
extra.plotmn = 1;
trunc.scheme = 'parall';
trunc.rad1 = 5;     
trunc.rad2 = 1;     

lattice  = create_lattice(a1,a2,128,128,trunc,extra); clearvars extra;

%% Declaring the birefringent geometric phase grating

epse = 1.7^2;
epso = 1.5^2;
epsm = diag([epse,epso]);
for iP = 1:size(lattice.P,1)
al = lattice.P(iP,1)/a1(1)*pi;
R = [cos(al),sin(al);-sin(al),cos(al)];
epsrot = R*epsm*R';
epsmask.xx(lattice.mnDir(iP,1),lattice.mnDir(iP,2)) = epsrot(1,1);
epsmask.xy(lattice.mnDir(iP,1),lattice.mnDir(iP,2)) = epsrot(1,2);
epsmask.yx(lattice.mnDir(iP,1),lattice.mnDir(iP,2)) = epsrot(2,1);
epsmask.yy(lattice.mnDir(iP,1),lattice.mnDir(iP,2)) = epsrot(2,2);
epsmask.z(lattice.mnDir(iP,1),lattice.mnDir(iP,2)) = epso;
end


plotpars.plotflag = 1;
plotpars.plotfignum = 50;
plotpars.N1t = 1; plotpars.N2t = 1;
plot_epsmask(lattice,epsmask,plotpars);           clearvars plotpars;

fteps = create_fteps(lattice,epsmask);

%% Declaring layer structure and solving RCWA

lambda = 4000;
dpg = lambda/2/(sqrt(epse)-sqrt(epso));
theta  = 0.01;
phi = 0;

%              superstr.    |  PG   |  subst.       
struc.d       = [10000         dpg     10000  ]; % 
struc.fts     = {             fteps          }; % material A 
struc.sigma   = [    0               0         ];    

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


%% Calculating the diffraction Jones matrices
extra.plotflag = 1;
extra.plotfignum = 21;

difdata = DJM_2d_basic_general(wm,lattice,struc,globpars,extra);

%% Comparing the result from RTA and from DJM
Ioutsub = 5; % using the Iout indexing
efficiency_RTA = RTAdata.T(  find( ~(RTAdata.IIsub - Ioutsub ) ) )
efficiency_DJM = sum( abs( difdata.Jtsupsub{1,Ioutsub} * [1;-1i]/sqrt(2) ).^2 )

%% Plotting the field

incpars.Iinc = 1;    
% Alternative expression to declare input channel in the (m,n) form.
%  difdata can be replaced with RTAdata, depending on what is available.
mninc = [0,0];  
switch inlayer
    case 'sup'  
    incpars.Iinc = find( ~(difdata.IIsup - find(~vecnorm((lattice.mn-mninc),2,2))) );
    case 'sub'
    incpars.Iinc = find( ~(difdata.IIsub - find(~vecnorm((lattice.mn-mninc),2,2))) );
end

incpars.inlayer = 'sup'; 
incpars.Esp     = [1,-1i];
fieldpars.type = 'a1_z';
fieldpars.N1 = 100;
fieldpars.N2 = 3000;
fieldpars.f = 0;

fieldpars.l = 2;
fieldpars.zz = 150;

fielddata = field_2d_basic_general(wm,lattice,struc,globpars,incpars,fieldpars);

figure(31)
h = pcolor(fielddata.X,fielddata.Z,real(fielddata.Ex)); hold on;
h.EdgeColor = 'none';
ax = gca; ax.YDir = 'reverse'; axis equal