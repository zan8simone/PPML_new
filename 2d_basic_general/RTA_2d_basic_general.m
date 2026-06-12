function RTAdata = RTA_2d_basic_general(wm,lattice,struc,globpars,incpars,extra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the reflection and transmission diffraction efficiencies 
% Calculates the absorption in the internal layers.
%
% 2-d structure, conical incidence
%
% Simone Zanotto, 2026
% -------------------------------------------------------------------------
%                                 INPUTS
%
% wm        -> the working matrices, as generated i.e. by   "solve_2d_basic_general.m"      [structure]
% lattice   -> the lattice data, as generated i.e. by "create_lattice.m"                    [structure]
% struc     -> the layer stack data  (see also "solve_2d_basic_general.m")                  [structure]
% globpars  -> the simulation global parameters (see also "solve_2d_basic_general.m")       [structure]
% extra     -> structure with fields
%                extra.plotflag    -> 1 or 0 to toggle plotting of open diffraction channel map
%                extra.plotfignum  -> the figure number where to plot
%
% incpars.inlayer  -> incidence layer ('sup' or 'sub')
% incpars.Esp      -> Jones vector of incident light in local sp basis.
%                       Not necessarily normalized.                      [2x1 complex]
%
%
% -------------------------------------------------------------------------
%                                 OUTPUTS
%
% RTAdata.Nsup  -> number of open diff. channels in superstrate                        [integer]
% RTAdata.IIsup -> a list of the I indices of open diffr. channels in the superstrate. [Nsup x 1 integer]
%
% RTAdata.Nsub  -> number of open diff. channels in suberstrate                        [integer]
% RTAdata.IIsub -> a list of the I indices of open diffr. channels in the suberstrate. [Nsub x 1 integer]
%
% RTAdata.abscosthsup ->  absolute value of the cosine of diffraction
%                             angles in superstrate (measured from z axis)             [Nsup x 1 real]
% RTAdata.abscosthsub ->  absolute value of the cosine of diffraction
%                             angles in suberstrate (measured from z axis)             [Nsub x 1 real]
%
% RTAdata.RR  -> intensity in reflected channels, indexed by the global index I    [lattice.N x 1 real]
% RTAdata.TT  -> intensity in transmitted channels, indexed by the global index I  [lattice.N x 1 real]
% RTAdata.AA  -> absorbed power in internal layers and conducting interfaces       [struc.L*2 + 1 real]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = size(wm.q,3)-2;
N = size(lattice.G,1);

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
RTAdata.IIsup = IIsup;

iprog = 1;
IIsub = [];
for I = 1:N
    if imag( wm.q(I,I,L+2) ) == 0 
        IIsub(iprog) = I;
        iprog = iprog+1;
    end
end
RTAdata.IIsub = IIsub;

Nsup = length(IIsup);
Nsub = length(IIsub);
RTAdata.Nsup = Nsup;
RTAdata.Nsub = Nsub;

Iinc = 1; % injecting light into the I = 1 channel, which is also the Isup = 1 or Isub = 1

%%%%%%% check if incidence channel is within bounds
switch incpars.inlayer
    case 'sup'
        if Iinc > Nsup
            error('I = 1, i.e. (m,n) = (0,0), is a closed channel. Check your setting of globpars.kparx and globpars.kpary')
        end
    case 'sub'
        if Iinc > Nsub
            error('I = 1, i.e. (m,n) = (0,0), is a closed channel. Check your setting of globpars.kparx and globpars.kpary')
        end
end


%%%%%%%%%%% plotting diffraction pattern
if extra.plotflag
figure(extra.plotfignum)

subplot(121) % diff pattern in superstrate
if imag(epssup) == 0
plot([-1,1]*k0*sqrt(epssup), [0,0],'k'); hold on;
plot([0,0],[-1,1]*k0*sqrt(epssup),'k'); hold on;
scatter(kparG(IIsup,1),kparG(IIsup,2),'ro'); hold on; axis equal
plot(cosd(0:5:360)*k0*sqrt(epssup),sind(0:5:360)*k0*sqrt(epssup),'k-')

textspac = min(k0*sqrt(epssup)/5,abs(lattice.g2(2))/5);
for i = 1:length(IIsup)
I = IIsup(i);
text(kparG(I,1),kparG(I,2)+textspac,...
    ['[',num2str(lattice.mn(I,1)),' ',num2str(lattice.mn(I,2)),']']);
text(kparG(I,1),kparG(I,2)+textspac*2,...
    ['I = ',num2str(I)]);
text(kparG(I,1),kparG(I,2)+textspac*3,...
    ['Isup = ',num2str(i)]);
end
else
text(-1,0,'Superstrate has Im(eps) != 0 \newline Diffraction orders are not defined')
axis([-1 1 -1 1])
end
title('Diffraction orders \newline in superstrate')
ax = gca; ax.YDir = 'reverse';


subplot(122) % diff pattern in substrate
if imag(epssub) == 0
plot([-1,1]*k0*sqrt(epssub), [0,0],'k'); hold on;
plot([0,0],[-1,1]*k0*sqrt(epssub),'k'); hold on;
scatter(kparG(IIsub,1),kparG(IIsub,2),'ro'); hold on; axis equal
plot(cosd(0:5:360)*k0*sqrt(epssub),sind(0:5:360)*k0*sqrt(epssub),'k-')

textspac = min(k0*sqrt(epssub)/5,abs(lattice.g2(2))/5);
for i = 1:length(IIsub)
I = IIsub(i);
text(kparG(I,1),kparG(I,2)+textspac,...
    ['[',num2str(lattice.mn(I,1)),' ',num2str(lattice.mn(I,2)),']']);
text(kparG(I,1),kparG(I,2)+textspac*2,...
    ['I = ',num2str(I)]);
text(kparG(I,1),kparG(I,2)+textspac*3,...
    ['Isub = ',num2str(i)]);
end

else
text(-1,0,'Substrate has Im(eps) != 0 \newline Diffraction orders are not defined')
axis([-1 1 -1 1])
end
title('Diffraction orders \newline in substrate')
ax = gca; ax.YDir = 'reverse';

end




%%%%%%%%%%%%%%%%%%%%%% local polarization basis & matrix Uin
Uintilde  = zeros(2*(Nsup+Nsub));

%%%%% superstrate
for i = 1:Nsup
I = IIsup(i);

%%% cosine of the diffraction angle
RTAdata.abscosthsup(i) =  abs( dot([kparG(I,:), wm.q(I,I,1)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,1)]) );

shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), wm.q(I,I,1)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uintilde((2*i-1):(2*i),(2*i-1):(2*i)) = U/( sqrt(sqrt(epssup))*sqrt( RTAdata.abscosthsup(i) ) );

end 


%%%%% substrate
for i = 1:Nsub
I = IIsub(i);

%%% cosine of the diffraction angle
RTAdata.abscosthsub(i) = abs( dot([kparG(I,:), wm.q(I,I,end)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,end)]) );

%%% input waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), -wm.q(I,I,end)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uintilde((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U/( sqrt(sqrt(epssub))*sqrt( RTAdata.abscosthsub(i) ) );

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
        Esp(2*Iinc-1 : 2*Iinc) = incpars.Esp;
    case 'sub'
        Esp(2*Nsup + 2*Iinc-1 : 2*Nsup + 2*Iinc) = incpars.Esp;
end


abinc = (Vin')*inv(Win)*Uintilde*Esp;



%%% a and b coefficients in all the layers
as = zeros(2*N,L+2); bs = zeros(2*N,L+2);
for l = 1:L+2
    as(:,l) = (eye(2*N) - wm.S2(:,:,l)*wm.S3(:,:,l))\(             wm.S1(:,:,l)*abinc(1:2*N) + wm.S2(:,:,l)*wm.S4(:,:,l)*abinc(2*N+1:end));
    bs(:,l) = (eye(2*N) - wm.S3(:,:,l)*wm.S2(:,:,l))\(wm.S3(:,:,l)*wm.S1(:,:,l)*abinc(1:2*N) +              wm.S4(:,:,l)*abinc(2*N+1:end));
end


%%% Incident wave flux
switch incpars.inlayer
    case 'sup'
    hpar = wm.phi(:,:,1)*abinc(1:2*N);
    epar =   wm.A(:,:,1)*abinc(1:2*N);
    case 'sub'
    hpar = wm.phi(:,:,end)*abinc(2*N+1:end);
    epar =  -wm.A(:,:,end)*abinc(2*N+1:end);
end

incflux = real(epar(N+1:end).*conj(hpar(N+1:end)) + epar(1:N).*conj(hpar(1:N))) ;


%%% Flux at the beginning of all layers (incl sup & sub)
fluxbeg = zeros(N,L+2);
for l = 1:L+2
hpar = wm.phi(:,:,l)*(as(:,l) + diag(exp(1i*diag(wm.q(:,:,l)*struc.d(l))))*bs(:,l));
epar =   wm.A(:,:,l)*(as(:,l) - diag(exp(1i*diag(wm.q(:,:,l)*struc.d(l))))*bs(:,l));

fluxbeg(:,l) = real(epar(N+1:end).*conj(hpar(N+1:end)) + epar(1:N).*conj(hpar(1:N))) ;
end


%%% Flux at the end of all layers (incl sup & sub)
fluxend = zeros(N,L+2);
for l = 1:L+2
hpar = wm.phi(:,:,l)*(diag(exp(1i*diag(wm.q(:,:,l)*struc.d(l))))*as(:,l) + bs(:,l));
epar =   wm.A(:,:,l)*(diag(exp(1i*diag(wm.q(:,:,l)*struc.d(l))))*as(:,l) - bs(:,l));

fluxend(:,l) = real(epar(N+1:end).*conj(hpar(N+1:end)) + epar(1:N).*conj(hpar(1:N))) ;
end

switch incpars.inlayer
    case 'sup'
        RR = (incflux-fluxend(:,1))/norm(incflux);
        TT = fluxbeg(:,end)/norm(incflux);
    case 'sub'
        RR = -(incflux-fluxbeg(:,end))/norm(incflux);
        TT = -fluxend(:,1)/norm(incflux);
end

AA = zeros(L+L+1,1);
for ll = 1:(L+L+1)
    if mod(ll,2) == 1
    AA(ll) = ( sum(fluxend(:,(ll-1)/2+1)) - sum(fluxbeg(:,(ll-1)/2+2)) )/norm(incflux);
    else
    AA(ll) = ( sum(fluxbeg(:,ll/2+1)) - sum(fluxend(:,ll/2+1)) )/norm(incflux);
    end
end

if abs(sum(AA)+sum(RR)+sum(TT)-1) > 1e-5 % Check energy conservation
    error('RTA:EnNotCons','Energy not conserved')
end

RTAdata.R = RR;
RTAdata.T = TT;
RTAdata.A = AA;
