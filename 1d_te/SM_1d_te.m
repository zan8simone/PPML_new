%%%
function [rl,rr,tlr,trl] = SM_1d_te(a,L,...
   epssup,epssub,epsyA,epsyB,sigma,f,d,...
   halfnpw,k0,kpar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the S-matrix.
%
% 1-d structure, TE polarization.
%
% Simone Zanotto, 2012-2026
% -------------------------------------------------------------------------
%                              INPUTS
%
% a              -> periodicity                                                 [real]
% L              -> number of layers excluded superstrate and substrate         [integer]
% epssup         -> superstrate permittivity       [real]               
% epssub         -> suberstrate permittivity       [complex]               
% epsyA          -> y- component of dielectric tensor for material A in the internal layers           [1xL, complex]
% epsyB          -> idem, for material B                                                              [1xL, complex]
% sigma          -> conductivity of interfaces                                       [1x(L+1), complex]
% f              -> duty cycle in the internal layers                                [1xL, real]
% d              -> thicknesses of the layers, including super- and sub-strate       [1x(L+2), real]
% halfnpw        -> half number of harmonic waves used in calculation   [integer]
% k0             -> wavevector in vacuum  (= 2*pi/lambda0)              [real]
% kpar           -> x-projection of the incident wavevector             [real]
%
% -------------------------------------------------------------------------
%                             OUTPUTS
%
%  The coefficients connect the amplitudes of the Ey fields, calculated at the interface between 
%   superstrate and first internal layer, or substrate and last internal layer.
%
% rl             -> reflection coefficient from superstrate to superstrate ("left reflectance")
% rr             -> reflection coefficient from substrate to substrate ("right reflectance")
% tlr            -> transmission coefficient from superstrate to substrate ("left to right")
% trl            -> transmission coefficient from substrate to superstrate ("right to left")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *the lengths (a,d) and inverse lengths (k0,kpar) must be set in the same
%  units (e.g. microns and inverse microns, respectively)
%
% *the units of sigma must be inverse Ohms.
%
% *has been tested only with real epssup and epssub
%
%
% *As well known, the RCWA has a single convergence parameter which is the truncation
%   order for the plane wave basis. Here, this is ruled by halfnpw, which is such that 
%   the complete basis is (-halfnpw,… 0, +halfnpw).
%    Always pay attention that the selected value for halfnpw guarantees a 
%    satisfactory convergence for the problem under analysis.
%
% *halfnpw = 0 sets the number of harmonic waves to 1, e.g. the scattering
%  matrix reduces to the ordinary 2x2 formalism for unpatterned multilayers
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = -halfnpw:halfnpw;
npw = size(n,2);
kx = (2*pi/a)*n + kpar;

% Initializing the work variables
q   = zeros(npw,npw,L+2);
phi = zeros(npw,npw,L+2);
A   = zeros(npw,npw,L+2);
						  
eps = zeros(npw,npw,L+2);


% Eigensolutions for super- and substrate
q(:,:,1)   = diag(sqrt_whittaker(epssup*k0^2 - kx.*kx));
phi(:,:,1) = eye(npw);
A(:,:,1)   = (eye(npw)*k0^2)/q(:,:,1);
							  
eps(:,:,1) = eye(npw)*epssup;

q(:,:,L+2)   = diag(sqrt_whittaker(epssub*k0^2 - kx.*kx));
phi(:,:,L+2) = eye(npw);
A(:,:,L+2)   = (eye(npw)*k0^2)/q(:,:,L+2);
								
eps(:,:,L+2) = eye(npw)*epssub;


% Eigensolutions for internal layers
if L > 0
for l = 1:L

   F = zeros(npw);
   
            for i = 1:npw % Fourier transform of square centered at the origin
            for j = 1:npw
                if (i == j)
                F(i,j) = f(l);
                else
                F(i,j) = sin(pi*f(l)*(n(i)-n(j)))/(pi*(n(i)-n(j)));
                end
            end
            end
 
   eps(:,:,l+1) = (  epsyB(l) -   epsyA(l))*F +   epsyA(l)*eye(npw);
   
   % Eigensolution calculation
   [phhi,qq] = eig(eps(:,:,l+1)*(k0^2*eye(npw)) - diag(kx)*diag(kx));
   q(:,:,l+1) = diag(sqrt_whittaker(diag(qq)));
   phi(:,:,l+1) = phhi;
   A(:,:,l+1) = (k0^2*eye(npw))*phhi/q(:,:,l+1);
   
end
end

% Forward propagation 
S1 = zeros(npw,npw,L+2); S2 = zeros(npw,npw,L+2);
S1(:,:,1) = eye(npw);
for l = 1:L+1
[S1(:,:,l+1),S2(:,:,l+1)] = smpropag_fw_cond(S1(:,:,l),S2(:,:,l),...
                                phi(:,:,l),phi(:,:,l+1),...
                                A(:,:,l),A(:,:,l+1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l+1)*d(l+1))),...
								sigma(l)*376.730/k0);
end

% Backward propagation
S3 = zeros(npw,npw,L+2); S4 = zeros(npw,npw,L+2);
S4(:,:,L+2) = eye(npw);
for ll = 1:L+1
    l = L+3-ll;
[S3(:,:,l-1),S4(:,:,l-1)] = smpropag_bw_cond(S3(:,:,l),S4(:,:,l),...
                                phi(:,:,l),phi(:,:,l-1),...
                                A(:,:,l),A(:,:,l-1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l-1)*d(l-1))),...
								sigma(l-1)*376.730/k0);
end

rl = -S3(halfnpw+1,halfnpw+1,1)   *exp(-1i*q(halfnpw+1,halfnpw+1,1)*d(1));        % left -> left reflectance
rr = -S2(halfnpw+1,halfnpw+1,L+2) *exp(-1i*q(halfnpw+1,halfnpw+1,L+2)*d(L+2));  % right -> right reflectance
tlr = S1(halfnpw+1,halfnpw+1,L+2) *exp(-1i*q(halfnpw+1,halfnpw+1,1)*d(1))    *A(halfnpw+1,halfnpw+1,L+2)/A(halfnpw+1,halfnpw+1,1);     % left -> right transmittance
trl = S4(halfnpw+1,halfnpw+1,1)   *exp(-1i*q(halfnpw+1,halfnpw+1,L+2)*d(L+2))*A(halfnpw+1,halfnpw+1,1)  /A(halfnpw+1,halfnpw+1,L+2);   % right -> left transmittance

