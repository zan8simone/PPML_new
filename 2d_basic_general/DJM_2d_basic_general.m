function difdata = DJM_2d_basic_general(wm,lattice,struc,globpars,extra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the Jones matrices for all pairs of open diffraction channels.
%
% 2-d structure, conical incidence
%
% Simone Zanotto, 2026
% -------------------------------------------------------------------------
%                              INPUTS
%
% wm        -> the working matrices, as generated i.e. by   "solve_2d_basic_general.m"      [structure]
% lattice   -> the lattice data, as generated i.e. by "create_lattice.m"                    [structure]
% struc     -> the layer stack data  (see also "solve_2d_basic_general.m")                  [structure]
% globpars  -> the simulation global parameters (see also "solve_2d_basic_general.m")       [structure]
% extra     -> structure with fields
%                extra.plotflag    -> 1 or 0 to toggle plotting of open diffraction channel map
%                extra.plotfignum  -> the figure number where to plot
%
% -------------------------------------------------------------------------
%                              OUTPUTS
%
% difdata.Nsup  -> number of open diff. channels in superstrate                        [integer]
% difdata.IIsup -> a list of the I indices of open diffr. channels in the superstrate. [Nsup x 1 integer]
%
% difdata.Nsub  -> number of open diff. channels in suberstrate                        [integer]
% difdata.IIsub -> a list of the I indices of open diffr. channels in the suberstrate. [Nsub x 1 integer]
%
% difdata.abscosthsup ->  absolute value of the cosine of diffraction
%                             angles in superstrate (measured from z axis)             [Nsup x 1 real]
% difdata.abscosthsub ->  absolute value of the cosine of diffraction
%                             angles in suberstrate (measured from z axis)             [Nsub x 1 real]
%
% difdata.Jsupsup -> Jones matrices of diffraction from super- to superstrate 
%                                         [structure containing Nsup x Nsup matrices, each 2x2 complex. 
%                                          The element {Isup,Jsup} describes diffraction 
%                                          from the Isup-th order to the Jsup-th order]
% difdata.Jsupsub -> Jones matrices of diffraction from super- to substrate
%                                         [structure containing Nsup x Nsub matrices, each 2x2 complex]
% difdata.Jsubsup -> Jones matrices of diffraction from sub- to superstrate
%                                         [structure containing Nsub x Nsup matrices, each 2x2 complex]
% difdata.Jsubsub -> Jones matrices of diffraction from sub- to substrate
%                                         [structure containing Nsub x Nsub matrices, each 2x2 complex]
%
% difdata.Jtsupsup (and analogous) -> rescaled Jones matrices (see manual)
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
difdata.IIsup = IIsup;

iprog = 1;
IIsub = [];
for I = 1:N
    if imag( wm.q(I,I,L+2) ) == 0 
        IIsub(iprog) = I;
        iprog = iprog+1;
    end
end
difdata.IIsub = IIsub;

Nsup = length(IIsup);
Nsub = length(IIsub);
difdata.Nsup = Nsup;
difdata.Nsub = Nsub;

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
for i = 1:Nsup
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
for i = 1:Nsub
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

%%%%%%%%%%%%%%%%%%%%%% local polarization basis & matrices Uin, Uout
Uin  = zeros(2*(Nsup+Nsub));
Uout = zeros(2*(Nsup+Nsub));

%%%%% superstrate
for i = 1:Nsup
I = IIsup(i);

%%% cosine of the diffraction angle
difdata.abscosthsup(i) =  abs( dot([kparG(I,:), wm.q(I,I,1)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,1)]) );
 
%%% input waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), wm.q(I,I,1)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uin((2*i-1):(2*i),(2*i-1):(2*i)) = U;
Uintilde((2*i-1):(2*i),(2*i-1):(2*i)) = U/( sqrt(sqrt(epssup))*sqrt( difdata.abscosthsup(i)) );

%%% output waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), -wm.q(I,I,1)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uout((2*i-1):(2*i),(2*i-1):(2*i)) = U;
Uouttilde((2*i-1):(2*i),(2*i-1):(2*i)) = U/( sqrt(sqrt(epssup))*sqrt( difdata.abscosthsup(i)) );

end 


%%%%% substrate
for i = 1:Nsub
I = IIsub(i);

%%% cosine of the diffraction angle
difdata.abscosthsub(i) = abs( dot([kparG(I,:), wm.q(I,I,end)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,end)]) );

%%% input waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), -wm.q(I,I,end)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uin((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U;
Uintilde((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U/( sqrt(sqrt(epssub))*sqrt( difdata.abscosthsub(i) ) );

%%% output waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), wm.q(I,I,end)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uout((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U;
Uouttilde((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U/( sqrt(sqrt(epssub))*sqrt( difdata.abscosthsub(i) ) );

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Wout matrix
Wout = zeros(2*Nsup+2*Nsub);
for Isub = 1:Nsub
I = IIsub(Isub);
Wout(2*Nsup + 2*Isub - 1, Isub       ) =  wm.A(N+I,  I,L+2); 
Wout(2*Nsup + 2*Isub - 1, Isub + Nsub) =  wm.A(N+I,N+I,L+2);
Wout(2*Nsup + 2*Isub    , Isub       ) = -wm.A(  I,  I,L+2);
Wout(2*Nsup + 2*Isub    , Isub + Nsub) = -wm.A(  I,N+I,L+2);
end
for Isup = 1:Nsup
I = IIsup(Isup);
Wout(2*Isup - 1, 2*Nsub + Isup       ) = -wm.A(N+I,  I,1) ;
Wout(2*Isup - 1, 2*Nsub + Isup + Nsup) = -wm.A(N+I,N+I,1) ;
Wout(2*Isup    , 2*Nsub + Isup       ) =  wm.A(  I,  I,1) ;
Wout(2*Isup    , 2*Nsub + Isup + Nsup) =  wm.A(  I,N+I,1) ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vin & Vout matrices
Vin = zeros(2*(Nsup+Nsub),4*N);
Vout = zeros(2*(Nsup+Nsub),4*N);
for Isup = 1:Nsup
I = IIsup(Isup);
Vin(Isup       ,     I ) = 1;
Vin(Isup + Nsup, N + I ) = 1;
Vout(2*Nsub + Isup       , 2*N + I ) = 1;
Vout(2*Nsub + Isup + Nsup, 3*N + I ) = 1;
end
for Isub = 1:Nsub
I = IIsub(Isub);
Vin(2*Nsup + Isub       , 2*N + I ) = 1;
Vin(2*Nsup + Isub + Nsub, 3*N + I ) = 1;
Vout(Isub       ,     I ) = 1;
Vout(Isub + Nsub, N + I ) = 1;
end


D   = inv(Uout)*Wout*Vout*[wm.S1(:,:,L+2), wm.S2(:,:,L+2); wm.S3(:,:,1), wm.S4(:,:,1)]*(Vin')*inv(Win)*Uin;
Dt  = inv(Uouttilde)*Wout*Vout*[wm.S1(:,:,L+2), wm.S2(:,:,L+2); wm.S3(:,:,1), wm.S4(:,:,1)]*(Vin')*inv(Win)*Uintilde;


difdata.Jsupsup = {}; difdata.Jsupsub = {}; difdata.Jsubsup = {}; difdata.Jsubsub = {};
difdata.Jtsupsup = {}; difdata.Jtsupsub = {}; difdata.Jtsubsup = {}; difdata.Jtsubsub = {};
% sup -> sup Jones matrices
for Isup = 1:Nsup
for Ksup = 1:Nsup
    difdata.Jsupsup{Isup,Ksup} = D( 2*Ksup-1 : 2*Ksup , 2*Isup-1 : 2*Isup );
    difdata.Jtsupsup{Isup,Ksup} = Dt( 2*Ksup-1 : 2*Ksup , 2*Isup-1 : 2*Isup );
end
end

% sup -> sub Jones matrices
for Isup = 1:Nsup
for Ksub = 1:Nsub
    difdata.Jsupsub{Isup,Ksub} = D(2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Isup-1 : 2*Isup );
    difdata.Jtsupsub{Isup,Ksub} = Dt(2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Isup-1 : 2*Isup );
end
end

% sub -> sup Jones matrices
for Isub = 1:Nsub
for Ksup = 1:Nsup
    difdata.Jsubsup{Isub,Ksup} = D( 2*Ksup-1 : 2*Ksup , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
    difdata.Jtsubsup{Isub,Ksup} = Dt( 2*Ksup-1 : 2*Ksup , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
end
end

% sub -> sub Jones matrices
for Isub = 1:Nsub
for Ksub = 1:Nsub
    difdata.Jsubsub{Isub,Ksub} = D( 2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
    difdata.Jtsubsub{Isub,Ksub} = Dt( 2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
end
end