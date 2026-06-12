function fielddata = field_2d_basic_general(wm,lattice,struc,globpars,incpars,fieldpars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the field
%
% 2-d structure, conical incidence
%
% Simone Zanotto, 2026
% incpars.Iinc     -> index of incidence channel (in [1 ... Nsup] or [1 ... Nsub])
% incpars.inlayer  -> incidence layer ('sup' or 'sub')
% incpars.Esp      -> Jones vector of incident light in local sp basis [2x1 complex]
%
% fieldpars.type   -> type of section ('a1_a2' [i.e., xy],   'a1_z',   'a2_z')
% fieldpars.N1     -> number of grid points along the first direction 
% fieldpars.N2     -> number of grid points along the second direction 
%
%  %% %%  %% %% parameters needed if type = 'a1_a2'
% fieldpars.l      -> layer number where the a1_a2-slice is calculated 
%                      (l=1-> superstrate, l = L+2 -> substrate)
% fieldpars.zz     -> z-coordinate within the layer at which the a1_a2-slice is calculated
%                      ( zz = 0 -> beginning of layer, zz = d(lz) -> end of layer)
%
%  %% %%  %% %% parameters needed if type = 'a1_z'
% fieldpars.f      -> fraction of basis vector a2 at which section is calculated
%                      ( in [0,1] )
%  %% %%  %% %% parameters needed if type = 'a2_z'
% fieldpars.f      -> fraction of basis vector a1 at which section is calculated
%                      ( in [0,1] )

% fielddata.X      -> matrix of X coordinates of grid
%   "    "  Y, Z   ->   "   "   Y, Z     "     "       "
% fielddata.Ex    -> matrix of Ex evaluated at grid coordinates
%   "    "   y, z      "   "   Ey, Ez   "     "       "
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = size(wm.q,3)-2;
N = lattice.N;

epssup = struc.epssup;
epssub = struc.epssub;

k0 = globpars.k0;
kparx = globpars.kparx;
kpary = globpars.kpary;

kx = lattice.G(:,1) + kparx;
ky = lattice.G(:,2) + kpary;

kparG = [kx,ky];


%%%%%%% determining the open diff orders into superstrate & substrate
iprog = 1;
IIsup = [];
for I = 1:N
    if imag( wm.q(I,I,1) ) == 0 
        IIsup(iprog) = I;
        iprog = iprog+1;
    end
end

iprog = 1;
IIsub = [];
for I = 1:N
    if imag( wm.q(I,I,L+2) ) == 0 
        IIsub(iprog) = I;
        iprog = iprog+1;
    end
end

Nsup = length(IIsup);
Nsub = length(IIsub);


%%%%%%% check if incidence channel is within bounds

switch incpars.inlayer
    case 'sup'
        if incpars.Iinc > Nsup
            error('Attempting incidence from a closed diffraction channel')
        end
    case 'sub'
        if incpars.Iinc > Nsub
            error('Attempting incidence from a closed diffraction channel')
        end
end



%%%%%%%%%%%%%%%%%%%%%% local polarization basis & matrix Uin
Uintilde  = zeros(2*(Nsup+Nsub));

%%%%% superstrate
for i = 1:Nsup
I = IIsup(i);

%%% cosine of the diffraction angle
abscosthsup =  abs( dot([kparG(I,:), wm.q(I,I,1)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,1)]) );

shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), wm.q(I,I,1)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uintilde((2*i-1):(2*i),(2*i-1):(2*i)) = U/( sqrt(sqrt(epssup))*sqrt( abscosthsup) );

end 


%%%%% substrate
for i = 1:Nsub
I = IIsub(i);

%%% cosine of the diffraction angle
abscosthsub = abs( dot([kparG(I,:), wm.q(I,I,end)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,end)]) );

%%% input waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), -wm.q(I,I,end)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uintilde((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U/( sqrt(sqrt(epssub))*sqrt( abscosthsub ) );

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Win matrix
Win  = zeros(2*Nsup+2*Nsub);
for Isup = 1:Nsup
I = IIsup(Isup);
Win(2*Isup - 1, Isup       ) =  wm.A(N+I,  I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));
Win(2*Isup - 1, Isup + Nsup) =  wm.A(N+I,N+I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));
Win(2*Isup    , Isup       ) = -wm.A(  I,  I,1) * exp(1i*wm.q(I,I,1)*struc.d(1)); 
Win(2*Isup    , Isup + Nsup) = -wm.A(  I,N+I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));
end
for Isub = 1:Nsub
I = IIsub(Isub);
Win(2*Nsup + 2*Isub - 1, 2*Nsup + Isub       ) = -wm.A(N+I,  I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
Win(2*Nsup + 2*Isub - 1, 2*Nsup + Isub + Nsub) = -wm.A(N+I,N+I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
Win(2*Nsup + 2*Isub    , 2*Nsup + Isub       ) =  wm.A(  I,  I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
Win(2*Nsup + 2*Isub    , 2*Nsup + Isub + Nsub) =  wm.A(  I,N+I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vin matrix
Vin = zeros(2*(Nsup+Nsub),4*N);
for Isup = 1:Nsup
I = IIsup(Isup);
Vin(Isup       ,     I ) = 1;
Vin(Isup + Nsup, N + I ) = 1;
end
for Isub = 1:Nsub
I = IIsub(Isub);
Vin(2*Nsup + Isub       , 2*N + I ) = 1;
Vin(2*Nsup + Isub + Nsub, 3*N + I ) = 1;
end


Esp = zeros(2*Nsup+2*Nsub,1);

switch incpars.inlayer
    case 'sup'
        Esp(2*incpars.Iinc-1 : 2*incpars.Iinc) = incpars.Esp;
    case 'sub'
        Esp(2*Nsup + 2*incpars.Iinc-1 : 2*Nsup + 2*incpars.Iinc) = incpars.Esp;
end


abinc = (Vin')*inv(Win)*Uintilde*Esp;

%%% a and b coefficients in all the layers
as = zeros(2*N,L+2); bs = zeros(2*N,L+2);
for l = 1:L+2
    as(:,l) = (eye(2*N) - wm.S2(:,:,l)*wm.S3(:,:,l))\(             wm.S1(:,:,l)*abinc(1:2*N) + wm.S2(:,:,l)*wm.S4(:,:,l)*abinc(2*N+1:end));
    bs(:,l) = (eye(2*N) - wm.S3(:,:,l)*wm.S2(:,:,l))\(wm.S3(:,:,l)*wm.S1(:,:,l)*abinc(1:2*N) +              wm.S4(:,:,l)*abinc(2*N+1:end));
end



%%% field calculation
switch fieldpars.type
    case 'a1_a2'

%%%%% begin direct lattice grid construction
for i1 = 1:fieldpars.N1
    for i2 = 1:fieldpars.N2
    P = lattice.a1*(i1-1)/fieldpars.N1 + lattice.a2*(i2-1)/fieldpars.N2;
    X(i1,i2) = P(1);
    Y(i1,i2) = P(2);
    end
end
clearvars i1 i2 P


l = fieldpars.l; zz = fieldpars.zz;
if zz > 0 && zz < struc.d(l) && (l-1) < (L+2)

epar =   wm.A(:,:,l)*(diag(exp(1i*diag(wm.q(:,:,l))* zz      ))*as(:,l) - ...
                      diag(exp(1i*diag(wm.q(:,:,l))*(struc.d(l)-zz)))*bs(:,l));
hpar =   wm.phi(:,:,l)*(diag(exp(1i*diag(wm.q(:,:,l))* zz      ))*as(:,l) + ...
                       diag(exp(1i*diag(wm.q(:,:,l))*(struc.d(l)-zz)))*bs(:,l));
ez   =   wm.etaz(:,:,l)*( diag(ky)*hpar(1:N) - diag(kx)*hpar(N+1:2*N) )/k0;

for i1 = 1:fieldpars.N1
for i2 = 1:fieldpars.N2
fielddata.Ex(i1,i2) =  ((exp(1i* (kx*X(i1,i2)+ky*Y(i1,i2)) )).')*epar(N+1:2*N) ;
fielddata.Ey(i1,i2) = -((exp(1i* (kx*X(i1,i2)+ky*Y(i1,i2)) )).')*epar(1:N) ;
fielddata.Ez(i1,i2) =  ((exp(1i* (kx*X(i1,i2)+ky*Y(i1,i2)) )).')*ez;
fielddata.Z(i1,i2) = sum(struc.d(1:l))-struc.d(l)+zz;
end
end

else
error('Invalid value for fieldpars.zz or fieldpars.l')
end

fielddata.X = X;
fielddata.Y = Y;

    case 'a1_z'

%%%%% begin direct lattice grid construction
z = linspace(0,sum(struc.d),fieldpars.N2);
for i1 = 1:fieldpars.N1
    for i2 = 1:fieldpars.N2
    P = lattice.a1*(i1-1)/fieldpars.N1  + lattice.a2*fieldpars.f;
    X(i1,i2) = P(1);
    Y(i1,i2) = P(2);
    fielddata.Z(i1,i2) = z(i2);
    end
end
clearvars i1 i2 P
Dsh = cumsum(struc.d);  % a shifted version of D 
D = [0,Dsh(1:end-1)];

for i2 = 1:fieldpars.N2
l = find(z(i2)<=Dsh,1,'first');
zz = z(i2)-D(l);

epar =   wm.A(:,:,l)*(diag(exp(1i*diag(wm.q(:,:,l))* zz      ))*as(:,l) - ...
                      diag(exp(1i*diag(wm.q(:,:,l))*(struc.d(l)-zz)))*bs(:,l));
hpar =   wm.phi(:,:,l)*(diag(exp(1i*diag(wm.q(:,:,l))* zz      ))*as(:,l) + ...
                        diag(exp(1i*diag(wm.q(:,:,l))*(struc.d(l)-zz)))*bs(:,l));
ez   =   wm.etaz(:,:,l)*( diag(ky)*hpar(1:N) - diag(kx)*hpar(N+1:2*N) )/k0;

for i1 = 1:fieldpars.N1
fielddata.Ex(i1,i2) =  ((exp(1i* (kx*X(i1,i2)+ky*Y(i1,i2)) )).')*epar(N+1:2*N) ;
fielddata.Ey(i1,i2) = -((exp(1i* (kx*X(i1,i2)+ky*Y(i1,i2)) )).')*epar(1:N) ;
fielddata.Ez(i1,i2) =  ((exp(1i* (kx*X(i1,i2)+ky*Y(i1,i2)) )).')*ez;
end
end
fielddata.X = X;
fielddata.Y = Y;
    

    case 'a2_z'
        error('case a2_z not yet implemented')
    otherwise
        error('invalid value for fieldpars.type')
end
