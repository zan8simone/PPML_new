function wm = solve_2d_basic_general(lattice,struc,globpars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the layer eigenvalue problems and the scattering matrix
% propagation
%
%  Mostly Whittaker-Culshaw (PRB 1999) & Liscidini-Sipe (PRB 2008) notation
%  No Fourier factorization rules
%  2d conical incidence
%    
% ========= INPUT
%
% lattice -> stucture as created by create_lattice(...)
%
% struc   ->   data about the layer structure. Said L the number of
%                               internal layers (ranging from 0 to infty):
%    struc.fts    -> a structure containing the Fourier transforms of
%                        permittivities (e.g. as created by create_fteps(...) ) [length L]
%    struc.d      -> vector of layer thicknesses        [L+2 real]
%    struc.sigma  -> vector of interface conductivities [L+1 complex]
%    struc.epssup -> superstrate permittivity [complex]
%    struc.epssub -> substrate permittivity [complex]
%
% globpars -> global simulation parameters
%    globpars.k0 -> the vacuum wavevector
%    globpars.kparx -> x-component of Bloch wavevector
%    globpars.kpary -> y-component of Bloch wavevector
%    
%  ======== OUTPUT
%
%  wm -> a structure containing the RCWA matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fts = struc.fts;
d   = struc.d;
sigma = struc.sigma;
epssup = struc.epssup;
epssub = struc.epssub;

k0 = globpars.k0;
kparx = globpars.kparx;
kpary = globpars.kpary;

if (length(d)-2) == length(fts) && (length(sigma)-1) == length(fts)
L = length(fts);
else
    error('Inconsistent length among fts, sigma, d');
end

npw = size(lattice.G,1);
kx = lattice.G(:,1) + kparx;
ky = lattice.G(:,2) + kpary;

% Initializing the work variables
wm.q   = zeros(2*npw,2*npw,L+2);
wm.phi = zeros(2*npw,2*npw,L+2);
wm.A   = zeros(2*npw,2*npw,L+2);

% Eigensolutions for super- and substrate
wm.q(1:npw,1:npw,1)     = diag(sqrt_whittaker(epssup*k0^2 - kx.*kx - ky.*ky));
wm.q(npw+1:2*npw,npw+1:2*npw,1) = wm.q(1:npw,1:npw,1);
wm.phi(:,:,1) = eye(2*npw);
eps = epssup;
wm.etaz(:,:,1) = eye(npw)/eps;

kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
wm.A(:,:,1)   = (1/k0)*((eye(2*npw)*k0^2 - kstorto)*wm.phi(:,:,1))/wm.q(:,:,1);

wm.q(1:npw,1:npw,L+2)     = diag(sqrt_whittaker(epssub*k0^2 - kx.*kx - ky.*ky));
wm.q(npw+1:2*npw,npw+1:2*npw,L+2) = wm.q(1:npw,1:npw,L+2);
wm.phi(:,:,L+2) = eye(2*npw);
eps = epssub;
wm.etaz(:,:,L+2) = eye(npw)/eps;

kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
wm.A(:,:,L+2)   = (1/k0)*((eye(2*npw)*k0^2 - kstorto)*wm.phi(:,:,L+2))/wm.q(:,:,L+2);


% Eigensolutions for internal layers
if L > 0
for l = 1:L

    fteps = fts{l};
    epsxx = fteps.xx;
    epsxy = fteps.xy;
    epsyx = fteps.yx;
    epsyy = fteps.yy;
    epsz  = fteps.z;   wm.etaz(:,:,l+1) = inv(epsz);


   % Eigensolution calculation
   kstorto = [(diag(ky)/epsz)*diag(ky), -(diag(ky)/epsz)*diag(kx);...
          -(diag(kx)/epsz)*diag(ky), (diag(kx)/epsz)*diag(kx)];
   kdritto = [diag(kx)*diag(kx), diag(kx)*diag(ky);...
              diag(ky)*diag(kx), diag(ky)*diag(ky)];
   epsilone = [epsyy, -epsyx; -epsxy, epsxx];

   [phhi,qq] = eig(epsilone*(eye(2*npw)*k0^2 - kstorto) ...
                    - kdritto);
   wm.q(:,:,l+1) = diag(sqrt_whittaker(diag(qq)));
   wm.phi(:,:,l+1) = phhi;
   wm.A(:,:,l+1) = (1/k0)*(eye(2*npw)*k0^2 - kstorto)*phhi/wm.q(:,:,l+1);
   
   
end
clearvars F
end



% Forward propagation 
wm.S1 = zeros(2*npw,2*npw,L+2); wm.S2 = zeros(2*npw,2*npw,L+2);
wm.S1(:,:,1) = eye(2*npw);
for l = 1:L+1
[wm.S1(:,:,l+1),wm.S2(:,:,l+1)] = smpropag_fw_cond(wm.S1(:,:,l),wm.S2(:,:,l),...
                                wm.phi(:,:,l),wm.phi(:,:,l+1),...
                                wm.A(:,:,l),wm.A(:,:,l+1),...
                                exp(1i*diag(wm.q(:,:,l)*d(l))),exp(1i*diag(wm.q(:,:,l+1)*d(l+1))),...
								sigma(l)*376.730);
end

% Backward propagation
wm.S3 = zeros(2*npw,2*npw,L+2); wm.S4 = zeros(2*npw,2*npw,L+2);
wm.S4(:,:,L+2) = eye(2*npw);
for ll = 1:L+1
    l = L+3-ll;
[wm.S3(:,:,l-1),wm.S4(:,:,l-1)] = smpropag_bw_cond(wm.S3(:,:,l),wm.S4(:,:,l),...
                                wm.phi(:,:,l),wm.phi(:,:,l-1),...
                                wm.A(:,:,l),wm.A(:,:,l-1),...
                                exp(1i*diag(wm.q(:,:,l)*d(l))),exp(1i*diag(wm.q(:,:,l-1)*d(l-1))),...
								sigma(l-1)*376.730);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

