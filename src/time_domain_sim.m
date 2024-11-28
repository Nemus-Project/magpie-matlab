function [outs,stresses,fs,fin] = time_domain_sim(rho,Evec,nux,Lvec,hvec,Nvec,KRmat,fmax,T,x_f,outMat,Nribs,beamParams,beamCoord,Nlump,Mlump,lumpCoord,dampVec,forceType,forceParams,LivePlot,RefreshRate,absPlot,FilmRec,cmap)


%-----------------------------------------------------------------------
%-- unpack constants
Lz      = Lvec(3) ;
hx      = hvec(1) ;
hy      = hvec(2) ;
Nx      = Nvec(1) ;
Ny      = Nvec(2) ;
Ex      = Evec(1) ;
Ey      = Evec(2) ;
Gxy     = Evec(3) ;
%-----------------------------------------------------------------------

nuy = nux*Ey/Ex ; % get second poisson ratio


%-----------------------------------------------------------------------
% set base stiffness and mass matrices
orthoBiHarm = orthoplate_biharm_build(Evec,nux,Lvec,hvec,Nvec,KRmat) ;
K = orthoBiHarm ;
M = speye((Nx+1)*(Ny+1))*rho*Lz   ;

% add ribs

for nR = 1 : Nribs

    Nbeam = beamParams(nR,5) ;
    x_coord_beam = linspace(beamCoord(nR,1),beamCoord(nR,2),Nbeam+1) ;
    y_coord_beam = linspace(beamCoord(nR+Nribs,1),beamCoord(nR+Nribs,2),Nbeam+1) ;
    RibDist = [] ;


    for nL = 1 : Nbeam + 1

        x_rib = [x_coord_beam(nL) ; y_coord_beam(nL)] ;
        RibDist = [RibDist,spreadinterp(x_rib,Nx,Ny,hx,hy,1)] ;

    end

    RibSpread = hx*hy*RibDist.' ;


    rhobeam = beamParams(nR,1) ;
    Ebeam = beamParams(nR,2) ;
    Ibeam = beamParams(nR,3) ;
    Abeam = beamParams(nR,4) ;
    hbeam = beamParams(nR,6) ;

    Bbeam = Ebeam*Ibeam*beam_biharmonic_build(Nbeam,hbeam*Nbeam,1) ;
    Mbeam = speye(Nbeam+1)*rhobeam*Abeam ;

    K = K+RibDist*Bbeam*RibSpread ;
    M = M+RibDist*Mbeam*RibSpread  ;

end

if Nlump
    % add lumped elements
    LumpDist = [] ;
    for nL = 1 : Nlump

        x_lump = [lumpCoord(nL) ; lumpCoord(nL+Nlump)] ;
        LumpDist = [LumpDist,spreadinterp(x_lump,Nx,Ny,hx,hy,1)] ;

    end

    LumpSpread = hx*hy*LumpDist.' ;

    M = M + LumpDist*diag(Mlump)*LumpSpread ;  M = sparse(M) ;
end

% quick trick to improve the matrices' numerical accuracy
M = 0.5*(M+M') ;
symK = 0.5*(K+K') ; 
skewK = 0.5*(K-K') ;
K = symK + skewK ;

[Q,Omsq] = eigs(K,M,(Nx+1)*(Ny+1),'smallestabs') ;
[~,indSort] = sort(diag((Omsq))) ;
Q = Q(:,indSort) ;
Om = sqrt(abs(diag(Omsq))) ;


% keep values up to fmax
fCur = -10 ;
indCur = 0 ;

while fCur < 0
        indCur = indCur + 1 ;
    fCur = Om(indCur) / 2 / pi ;
end
Om = Om(indCur + 1:end) ;
Q = Q(:,indCur + 1:end) ;

fCur = -10 ;
indCur = 0 ;

while fCur < fmax 
    indCur = indCur + 1 ;
    fCur = Om(indCur) / 2 / pi ;
end
Om = Om(1:indCur) ;
Q = Q(:,1:indCur) ;


Nmodes = indCur ;

%-----------------------------------------------------------------------



%-----------------------------------------------------------------------
%-- time-domain parameters
k = 2/Om(end) ;
fs = round(1/k) ; % sample rate
k = 1 / fs ;      % time step

Ts = round(fs*T) ; % number of simulation steps
tv = (0:Ts-1)*k ;   % time vector
fv = (0:Ts-1)*fs/Ts ; % frequency vector
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-- forcing parameters
if forceType == 1
    tw = forceParams(1) ;  % total contact time [s]
    F0 = forceParams(2) ;  % largest forcing amp [N]
    b = forceParams(3) ;   % noise modulation

    noise = rand(1,Ts) ;
    tws = floor(tw*fs) ;

    % raised cosine
    rc  = zeros(1,Ts) ;
    rc(1:tws) = 0.5*F0*(1 - cos(2*pi*(0:tws-1)/tws)) ;
    fin   = rc.*(1+b*noise) ;
elseif forceType == 2
    forceFreq = forceParams(1) ;  % total contact time [s]
    F0 = forceParams(2) ;  % largest forcing amp [N]
    b = forceParams(3) ;   % noise modulation

    noise = rand(1,Ts) ;

    % raised cosine
    fin = F0*sin(2*pi*forceFreq*(0:Ts-1)*k);
    fin   = fin.*(1+b*noise) ;
else
    [fin,fsIn] = audioread(forceParams) ;
    fin = resample(fin,fs,fsIn) ;
    if Ts > length(fin) 
    fin = [fin; zeros(Ts-length(fin),1)] ;
    end
end
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-- losses
OmDamp1 = 0 ;
OmDamp2 = 2*pi*dampVec(3) ;
t60_1   = dampVec(1) ;
t60_2   = dampVec(2) ;
dOmSq   = (OmDamp2^2 - OmDamp1^2) ;
alpha   = 3*log(10)/dOmSq * (OmDamp2^2/t60_1 - OmDamp1^2/t60_2) ; % rayleigh coefficient 1
beta    = 3*log(10)/dOmSq * (1/t60_2 - 1/t60_1) ;  % rayleigh coefficient 2
sig     = alpha + beta*Om.^2 ; % loss factor
%sig     = [sig;CLump] ;
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-- I/O init

%input
Jin   = spreadinterp(x_f,Nx,Ny,hx,hy,1) ;
% outputs
Nouts = length(outMat(:,1)) ;
Jouts = [] ;
for no = 1 : Nouts

    o_p = outMat(no,:) ;

    Jouts = [Jouts; [spreadinterp(o_p,Nx,Ny,hx,hy,2)]] ;
end
%-----------------------------------------------------------------------


%-----------------------------------------------------------------------
%-- matrix init
Dxx = DxxBuild(Nx,Ny,hx,'blk') ;
Dyy = DyyBuild(Nx,Ny,hy,'blk') ;
Dxy = DxyBuild(Nx,Ny,hx,hy,'blk') ;

interpouts = Jouts*Q ;  % displacement interpolator
interpstrainx = -0.5*Lz*Jouts*Dxx*Q ;  % strain x interpolator
interpstrainy = -0.5*Lz*Jouts*Dyy*Q ;  % strain y interpolator
interpstrainxy = -Lz*Jouts*Dxy*Q ;  % strain xy interpolator

Jin  = Q \ (M \ Jin) ; % input dirac
%-----------------------------------------------------------------------


%-----------------------------------------------------------------------
%-- state-vector init (modal coordinates)
vm      = zeros(Nmodes,1) ;
v0      = zeros(Nmodes,1) ;
P0      = 2*cos(Om*k).*exp(-sig*k) ;
Pm      = exp(-2*sig*k) ;
F       = k^2*exp(-sig*k).*Jin ;
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%-- memory init
w = zeros(Nouts,Ts) ;
strainx = zeros(Nouts,Ts) ;
strainy = zeros(Nouts,Ts) ;
strainxy = zeros(Nouts,Ts) ;
%-----------------------------------------------------------------------


%-----------------------------------------------------------------------
%-- plot init
close all ;

xax = (0 : Nx)*hx ;
yax = (0 : Ny)*hy ;

[xax,yax] = meshgrid(xax,yax) ;  % axis for 3D plots
colormap(cmap) ;
%-----------------------------------------------------------------------


%-----------------------------------------------------------------------
%-- main loop
if LivePlot == 1
    if FilmRec == 1
        vidObj = VideoWriter('DynamicPlots');
        vidObj.FrameRate = 10;
        open(vidObj) ;
    end
end

tic
for n = 1 : Ts

    vp = P0.*v0 - Pm.*vm + F*fin(n) ; % modal scheme
    w(:,n) = interpouts*v0 ; % displacement output
    % strains
    strainx(:,n)  = interpstrainx*v0 ;
    strainy(:,n)  = interpstrainy*v0 ;
    strainxy(:,n) = interpstrainxy*v0 ;

    if LivePlot == 1 %3D plot of stress fields and disp,vel,acc

        if n == 1 || mod(n,RefreshRate) == 0

            disp = Q*vp ;
            vel  = Q*(vp-vm)/2/k ;
            acc  = Q*(vp-2*v0+vm)/k^2 ;

            ex  = -0.5*Lz*Dxx*disp ;
            ey  = -0.5*Lz*Dyy*disp ;
            gm  = -Lz*Dxy*disp ;

            sigx = Ex/(1-nux*nuy)*(ex+nuy*ey) ;
            sigy = Ey/(1-nux*nuy)*(ey+nux*ex) ;
            tau  = Gxy*gm ;

            figure; colormap(cmap) ;


            subplot(2,3,1)
            animated_plot(sigx,1e6,xax,yax,Nvec,Lvec,absPlot,'$\sigma_x$ (Mpa)',Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f) ;

            subplot(2,3,2)
            animated_plot(sigy,1e6,xax,yax,Nvec,Lvec,absPlot,'$\sigma_y$ (Mpa)',Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f) ;

            subplot(2,3,3)
            animated_plot(tau,1e6,xax,yax,Nvec,Lvec,absPlot,'$\tau$ (Mpa)',Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f) ;

            subplot(2,3,4)
            animated_plot(disp,1e-3,xax,yax,Nvec,Lvec,absPlot,'disp (mm)',Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f) ;

            subplot(2,3,5)
            animated_plot(vel,1,xax,yax,Nvec,Lvec,absPlot,'vel (m/s)',Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f) ;

            subplot(2,3,6)
            animated_plot(acc,1e3,xax,yax,Nvec,Lvec,absPlot,'acc (km/s)',Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f) ;

            tit = sprintf('Time = %1.1f ms', k*(n-1)/1e-3) ;
            sgtitle(tit,'interpreter','latex') ;

            drawnow ;

            if FilmRec == 1
                frame = getframe(gcf);
                writeVideo(vidObj,frame)
            end

        end

    end

    vm = v0 ; v0 = vp ; % update

end
toc

if LivePlot == 1
    if FilmRec == 1
        close(vidObj) ;
    end
end
%- get outputs

vels = diff(w,1,2)/k ;
vels = [vels,vels(:,end)] ;
accs = diff(vels,1,2)/k ;
accs = [accs,accs(:,end)] ;

outs = [w;vels;accs] ;

sigx = Ex/(1-nux*nuy)*(strainx+nuy*strainy) ;
sigy = Ey/(1-nux*nuy)*(strainy+nux*strainx) ;
tau  = Gxy*strainxy ;

stresses = [sigx; sigy; tau] ;

% plot outputs
close all

figure()
subplot(3,1,1)
plot(tv,w/1e-3) ;
xlabel('$t$ (s)', 'interpreter', 'latex') ;
ylabel('$u$ (mm)', 'interpreter', 'latex') ;
title('displacement','interpreter', 'latex')

subplot(3,1,2)
plot(tv,vels) ;
xlabel('$t$ (s)', 'interpreter', 'latex') ;
ylabel('$\dot u$ (m/s)', 'interpreter', 'latex') ;
title('velocity','interpreter', 'latex')

subplot(3,1,3)
plot(tv,accs) ;
xlabel('$t$ (s)', 'interpreter', 'latex') ;
ylabel('$\ddot u$ (m/s$^2$)', 'interpreter', 'latex') ;
title('acceleration','interpreter', 'latex')


figure()
subplot(3,1,1)
plot(tv,sigx/1e6) ;
xlabel('$t$ (s)', 'interpreter', 'latex') ;
ylabel('$\sigma_x$ (MPa)', 'interpreter', 'latex') ;
title('normal stress $x$','interpreter', 'latex')

subplot(3,1,2)
plot(tv,sigy/1e6) ;
xlabel('$t$ (s)', 'interpreter', 'latex') ;
ylabel('$\sigma_y$ (MPa)', 'interpreter', 'latex') ;
title('normal stress $y$','interpreter', 'latex')

subplot(3,1,3)
plot(tv,tau/1e6) ;
xlabel('$t$ (s)', 'interpreter', 'latex') ;
ylabel('$\tau$ (MPa)', 'interpreter', 'latex') ;
title('shear stress','interpreter', 'latex')


figure()
plot(fv,20*log10(abs(fft(vels')))) ; xlim([0,fs/2]) ;
xlabel('$f$ (Hz)', 'interpreter', 'latex') ;
ylabel('$|\hat{\dot u}|$ (dB)', 'interpreter', 'latex') ;
title('velocity spectra','interpreter', 'latex')




end