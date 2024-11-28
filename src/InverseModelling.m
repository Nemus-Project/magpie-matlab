clear all
close all


%------------- Training Phase: Obtain the linear modal paramters a b c

%-------------------------------------------------------------------------
% CUSTOM PARAMETERS (user can change these)


%--------------------
%-- plate parameters
rho      = 473.9 ;
Lx       = 0.223 ;
Ly       = 0.114 ;
Lz       = 0.003 ;

%-- guessed values 
Ex0      = 10.7e9 ;   
Ey0      = 716e6 ;
Gxy0     = 500e6 ;
nux0     = 0.51 ;

%-- elastic constants around the edges
KRmat = [0e1,0e13; %Kx0 Rx0
    1e13, 1e13; % K0y R0y
    0e13, 0e13; % KxL RxL
    0e13, 0e13] ; % KLy RLy

%-- measured experimental frequencies (Hz)
ExpFreqs = [52
    98
    311
    337
    398
    637] ;
%--------------------


%--------------------
%-- general parameters
fmax = 2000 ; %-- maximum frequency to be computed
ppw  = 100 ; % points per wavelength at maximum frequency. Choose 3 <= ppw
Nmeas = 6 ; % total number of modes required
%--------------------

%--------------------
%-- braces parameters
Nribs = 0 ; % number of braces
Eb = [].' ; % youngs moduli
Lzb = [].' ; % thicknesses
bb = [].' ; % width cross section
rhob = [].' ; % densities

% rib coordinates along x (start and end) AS A FRACTION OF Lx
x_beam_coord = ...
    [] ;

% rib coordinates along y (start and end) AS A FRACTION OF Ly
y_beam_coord = ...
    [] ;

%--------------------
%-- static loads and stiffeners parameters

Nlump = 0 ; % number of lumped elements
x_lump_coord = [].' ; % x coordinates of lumped elements
y_lump_coord = [].' ; % y coordinates of lumped elements
KLump  = [].' ;
MLump  = [].' ;
%--------------------

%--------------------
%-- plot parameters parameters
cmap = cmaps(1) ; % select colormap 1 = RedBlue, 2 = GreenPurple, 3 = OrangeGreen, 4 = PurpleOrange
absPlot = 0 ;
NN = 1 ; % mode number to plot (plot will display modes between NN and NN + 8)
%--------------------

%--------------------
%-- MAC parameters
MAC_threshold = 0.99 ;
%--------------------

% END CUSTOM PARAMETERS
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% DERIVED PARAMETERS (user cannot change these)

%-- plate grid spacings and grid points
Evec = [Ex0,Ey0,Gxy0] ;
Lvec = [Lx,Ly,Lz] ;
facx = sqrt(2*pi)*sqrt(sqrt(Evec(1)*Lvec(3)^2/12/rho)) ;
facy = sqrt(2*pi)*sqrt(sqrt(Evec(2)*Lvec(3)^2/12/rho)) ;
hx = facx/ppw/sqrt(fmax) ; % this applies the formula hx*ppw = cx/fmax (\lambda f = c)
hy = facy/ppw/sqrt(fmax) ; % this applies the formula hy*ppw = cy/fmax (\lambda f = c)
Nx   = round(Lvec(1)/hx) ; Ny = round(Lvec(2)/hy) ;
Nvec = [Nx,Ny] ; % grid points
hvec = [Lvec(1)/Nvec(1), Lvec(2)/Nvec(2)] ; % grid spacings


%-- braces parameters
Ib = bb.*Lzb.^3/12 ; % moments of inertia
Ab = bb.*Lzb ; %  areas
hb = sqrt(sqrt(Eb.*Ib./rhob./Ab))*sqrt(2*pi/fmax)/ppw ; %  grid spacings in terms of ppw
Lb = [] ; Nb = [] ;

if Nribs
    Lb = sqrt((x_beam_coord(:,1)-x_beam_coord(:,2)).^2 + (y_beam_coord(:,1)-y_beam_coord(:,2)).^2) ;
    Nb = round(Lb./hb) ; % number of points
    hb = Lb./Nb ; % actual grid spacings
end

beamParams = [rhob,Eb,Ib,Ab,Nb,hb] ;
beamCoord = [x_beam_coord; y_beam_coord] ;


%-- lumped elements parameters
lumpParams = [KLump,MLump] ;
lumpCoord = [x_lump_coord; y_lump_coord] ;


%-- plot parameters
close all
xax = linspace(0,Lx,Nx+1) ;
yax = linspace(0,Ly,Ny+1) ;
[X,Y] = meshgrid(xax,yax) ;
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% EIGENVALUE PROBLEM
[Om,Q] = freq_domain_sim(rho,Evec,nux0,Lvec,hvec,Nvec,KRmat,fmax,Nmeas,Nribs,beamParams,beamCoord,Nlump,lumpParams,lumpCoord) ;
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% PLOTS
modal_plotter(cmap,6,NN,X,Y,Q,Lvec,Nvec,Nribs,Nlump,beamParams,beamCoord,lumpCoord)
%-------------------------------------------------------------------------

QRef = Q ;
OmRef = Om ;



scaleVec = [0.8,0.9,1,1.1,1.2] ;
scaleMat = zeros(2,5^2) ;
ind = 0 ;
indRemove = [] ;
OmMat = zeros(6,25) ;
for p = 1 : 5
    for q = 1 : 5
        ind = ind + 1
        scaleMat(:,ind) = [scaleVec(p);scaleVec(q)] ;
        Ex      = Ex0*scaleVec(p) ;
        Ey      = Ey0*scaleVec(q) ;
        Evec = [Ex,Ey,Gxy0] ;
        [Om,Q] = freq_domain_sim(rho,Evec,nux0,Lvec,hvec,Nvec,KRmat,fmax,2*Nmeas,Nribs,beamParams,beamCoord,Nlump,lumpParams,lumpCoord) ;
        Qtemp = [] ;
        MAC = [] ;
        for pMAC = 1 : Nmeas
            modeRef = QRef(:,pMAC) ;
            maxMAC = -10 ; MAC_index = 0 ;
            for qMAC = 1 : 2*Nmeas

                modeCur = Q(:,qMAC) ;

                MACcur =  sqrt(abs(modeRef.' * modeCur)^2 / (modeRef.' * modeRef) / (modeCur.' * modeCur)) ;

                if MACcur > maxMAC
                    maxMAC = MACcur ;
                    MAC_index = qMAC ;
                end
            end

            OmMat(pMAC,ind) = Om(MAC_index) ;
            Qtemp = [Qtemp,Q(:,MAC_index)] ;
            MAC = [MAC,maxMAC] ;

        end

        if mean(MAC) < MAC_threshold
            indRemove = [indRemove,ind] ;
        end
        %-------------------------------------------------------------------------
        % PLOTS
        close all
        modal_plotter(cmap,6,NN,X,Y,Q,Lvec,Nvec,Nribs,Nlump,beamParams,beamCoord,lumpCoord)
        sttit = sprintf('Test Index = %d, Avg. MAC = %1.2f',ind,mean(MAC)) ;
        sgtitle(sttit,'interpreter','latex') ;
        drawnow ; pause(0.1) ;
        %-------------------------------------------------------------------------

    end
end

OmMat(:,[indRemove]) = [] ;
scaleMat(:,[indRemove]) = [] ;

Nsets = length(OmMat(1,:)) ;

%-- init
DxVec = zeros(Nsets,1) ;
DyVec = zeros(Nsets,1) ;
DsVec = zeros(Nsets,1) ;
om1sq = zeros(Nsets,1) ;
om2sq = zeros(Nsets,1) ;
om3sq = zeros(Nsets,1) ;
om4sq = zeros(Nsets,1) ;
om5sq = zeros(Nsets,1) ;
om6sq = zeros(Nsets,1) ;
rxVec = zeros(Nsets,1) ;
ryVec = zeros(Nsets,1) ;

for n = 1 : Nsets


    Ex      = Ex0*scaleMat(1,n) ;
    Ey      = Ey0*scaleMat(2,n) ;
    Gxy     = Gxy0 ;

    nux     = nux0 ;
    nuy     = Ey/Ex*nux ;

    Dx            = Ex*Lz^3/12/(1-nux*nuy) ;
    Dy            = Ey*Lz^3/12/(1-nux*nuy) ;
    Ds            = Gxy*Lz^3/3 ;

    DxVec(n)      = Dx ;
    DyVec(n)      = Dy ;
    DsVec(n)      = Ds ;
    rxVec(n)      = Dx/Ds ;
    ryVec(n)      = Dy/Ds ;
    om1sq(n)      = OmMat(1,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om2sq(n)      = OmMat(2,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om3sq(n)      = OmMat(3,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om4sq(n)      = OmMat(4,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om5sq(n)      = OmMat(5,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;
    om6sq(n)      = OmMat(6,n).^2/Ds*rho*Lz*Lx^2*Ly^2 ;

end



omsq = [om1sq,om2sq,om3sq,om4sq,om5sq,om6sq] ;



FitMat = zeros(4,6) ; %-- this collects the best fit paramters for the six modes
%-- first row: a ; second row: b ; third row: c.

for n = 1 : 6

    linfit = fit([rxVec,ryVec],omsq(:,n),'poly11','Lower', [0 0 0]) ;
    %linfit = fit([rxVec,ryVec],omsq(:,n),'poly11') ;

    aCur = linfit.p10 ;
    bCur = linfit.p01 ;
    cCur = linfit.p00 ;


    FitMat(1,n) = aCur ;
    FitMat(2,n) = bCur ;
    FitMat(3,n) = cCur ;

    %-- compute Rsquared
    SStot = 0 ; SSres = 0 ; meany = mean(omsq(:,n)) ;

    for nFit = 1 : Nsets
        rx        = rxVec(nFit) ;
        ry        = ryVec(nFit) ;
        yfit      = aCur*rx + bCur*ry + cCur ;
        y         = omsq(nFit,n) ;
        SStot     = SStot + (y-meany).^2;                    % Total Sum-Of-Squares
        SSres     = SSres + (y-yfit).^2 ;                    % Residual Sum-Of-Squares

    end
    Rsq       = 1-SSres/SStot  ;
    FitMat(4,n) = Rsq ;

end


%------------- End Training Phase ........................................

%--------------- Plate: Measured spruce plate

indCell = {
    [1 2 3]
    [1 2 4]
    [1 2 5 ]
    [1 2 6 ]
    [1 3 4 ]
    [1 3 5 ]
    [1 3 6 ]
    [1 4 5 ]
    [1 4 6 ]
    [1 5 6 ]
    [2 3 4]
    [2 3 5]
    [2 3 6]
    [2 4 5]
    [2 4 6]
    [2 5 6]
    [3 4 5]
    [3 4 6]
    [3 5 6]
    [4 5 6]
    [1 2 3 4]
    [1 2 3 5]
    [1 2 3 6]
    [1 2 4 5]
    [1 2 4 6]
    [1 2 5 6]
    [1 3 4 5]
    [1 3 4 6]
    [1 3 5 6]
    [1 4 5 6]
    [2 3 4 5]
    [2 3 4 6]
    [2 3 5 6]
    [2 4 5 6]
    [3 4 5 6]
    [1 2 3 4 5]
    [1 2 3 4 6]
    [1 2 3 5 6]
    [1 2 4 5 6]
    [1 3 4 5 6]
    [2 3 4 5 6]
    [1 2 3 4 5 6]} ;
Ntot = 42 ;

errDx = zeros(Ntot,1) ;
errDy = zeros(Ntot,1) ;
errDs = zeros(Ntot,1) ;
DxLS  = zeros(Ntot,1) ;
DyLS  = zeros(Ntot,1) ;
DsLS  = zeros(Ntot,1) ;
condNumb = zeros(Ntot,1) ;

rho      = 473.9 ;
Lx       = 0.223 ;
Ly       = 0.114 ;
Lz       = 0.003 ;
nux      = nux0 ;
nuy      = Ey0/Ex0*nux ;


Om0=2*pi*ExpFreqs;

for nCase = 1 : Ntot


    Nselect = indCell(nCase,:) ;
    Nselect = cell2mat(Nselect) ;


    a        = (FitMat(1,Nselect)).' ;
    b        = (FitMat(2,Nselect)).' ;
    c        = (FitMat(3,Nselect)).' ;
    rsquared = (FitMat(4,Nselect)).' ;



    Om = Om0(Nselect) ;


    psi  = Om.^2 ;


    eta = (rho*Lz*Lx^2*Ly^2)^(-1) ;

    X = [eta*a eta*b eta*c] ;
    LSmat = (X.' * X) \ (X.') ;
    condNumb(nCase) = cond(LSmat) ;


    temp = (X.' * X) \ (X.'  * psi) ; % least-square solution


    DxLS(nCase)     = temp(1) ;
    DyLS(nCase)     = temp(2) ;
    DsLS(nCase)     = temp(3) ;


end


disp('Test Experimental --------------------------------')


%Obtained Ex, Ey, Gxy



digits(3) ;
disp('Ex Ey Gxy')

ExLSall           = (DxLS*(12*(1-nux*nuy)))/Lz^3 ;
EyLSall           = (DyLS*(12*(1-nux*nuy)))/Lz^3 ;
GxyLSall          = 3/Lz^3*(DsLS) ;

ExLS              = [(1:Ntot)',ExLSall] ;
EyLS              = [(1:Ntot)',EyLSall] ;
GxyLS             = [(1:Ntot)',GxyLSall] ;

EMat   = [ExLSall EyLSall GxyLSall] ;

ExLSNeg              = find(ExLS(:,2)<0) ;   ExLS(ExLSNeg,:) = [] ;
EyLSNeg              = find(EyLS(:,2)<0) ;   EyLS(EyLSNeg,:) = [] ;
GxyLSNeg             = find(GxyLS(:,2)<0) ;  GxyLS(GxyLSNeg,:) = [] ;

gigiEx   = isoutlier(ExLS(:,2),"quartiles") ;
gigiEy   = isoutlier(EyLS(:,2),"quartiles") ;
gigiGxy  = isoutlier(GxyLS(:,2),"quartiles") ;

rEx      = find(gigiEx) ;
indOutEx = ExLS(rEx,1) ;
ExLS(rEx,:)  = [] ;

rEy      = find(gigiEy) ;
indOutEy = EyLS(rEy,1) ;
EyLS(rEy,:)  = [] ;

rGxy      = find(gigiGxy) ;
indOutGxy = GxyLS(rGxy,1) ;
GxyLS(rGxy,:)  = [] ;

meanEx = mean(ExLS(:,2)) ; stdEx = std(ExLS(:,2)) ;
meanEy = mean(EyLS(:,2)) ; stdEy = std(EyLS(:,2)) ;
meanGxy = mean(GxyLS(:,2)) ; stdGxy = std(GxyLS(:,2)) ;

lengthVec = [length(ExLS), length(EyLS), length(GxyLS)] ;
Results_ElasticConstants   = vpa([meanEx meanEy meanGxy])
Standard_Deviations_Percent    = vpa([stdEx/meanEx*100 stdEy/meanEy*100 stdGxy/meanGxy*100])


Nstd = 8 ;

subplot(1,3,1)
plot((0:Ntot+1),ones(Ntot+2,1)*meanEx/1e9,'--k') ; xlim([0,Ntot+1]); ylim([meanEx-Nstd*stdEx,meanEx+Nstd*stdEx]/1e9);
x = [0 Ntot+1 Ntot+1 0];
y = [meanEx-stdEx meanEx-stdEx meanEx+stdEx meanEx+stdEx]/1e9;
patch(x,y,'red','FaceAlpha',0.2,'LineStyle','none'); hold on ;
plot(ExLS(:,1),ExLS(:,2)/1e9,'linestyle','none','marker','o','color','k');
%plot(ExLSNeg,ExLSall(ExLSNeg),'linestyle','none','marker','*','color','k') ;
plot(indOutEx,ExLSall(indOutEx)/1e9,'linestyle','none','marker','+','color','k');
legend('mean','std','kept','outlier','Location','southeast');
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
ylabel('$\tilde E_x$ (GPa)','interpreter','latex') ;
xlabel('mode set','interpreter','latex') ;

subplot(1,3,2)
plot((0:Ntot+1),ones(Ntot+2,1)*meanEy/1e6,'--k') ; xlim([0,Ntot+1]); ylim([meanEy-Nstd*stdEy,meanEy+Nstd*stdEy]/1e6)
x = [0 Ntot+1 Ntot+1 0];
y = [meanEy-stdEy meanEy-stdEy meanEy+stdEy meanEy+stdEy]/1e6;
patch(x,y,'blue','FaceAlpha',0.2,'LineStyle','none'); hold on ;
plot(EyLS(:,1),EyLS(:,2)/1e6,'linestyle','none','marker','o','color','k');
%plot(EyLSNeg,EyLSall(EyLSNeg),'linestyle','none','marker','*','color','k') ;
plot(indOutEy,EyLSall(indOutEy)/1e6,'linestyle','none','marker','+','color','k');
legend('mean','std','kept','outlier','Location','southeast');
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
ylabel('$\tilde E_y$ (MPa)','interpreter','latex') ;
xlabel('mode set','interpreter','latex') ;


subplot(1,3,3)
plot((0:Ntot+1),ones(Ntot+2,1)*meanGxy/1e6,'--k') ; xlim([0,Ntot+1]); ylim([meanGxy-Nstd*stdGxy,meanGxy+Nstd*stdGxy]/1e6)
x = [0 Ntot+1 Ntot+1 0];
y = [meanGxy-stdGxy meanGxy-stdGxy meanGxy+stdGxy meanGxy+stdGxy]/1e6;
patch(x,y,'green','FaceAlpha',0.2,'LineStyle','none'); hold on ;
plot(GxyLS(:,1),GxyLS(:,2)/1e6,'linestyle','none','marker','o','color','k');
%plot(GxyLSNeg,GxyLSall(GxyLSNeg),'linestyle','none','marker','*','color','k') ;
plot(indOutGxy,GxyLSall(indOutGxy)/1e6,'linestyle','none','marker','+','color','k');
legend('mean','std','kept','outlier','Location','southeast');
set(gca,'TickLabelInterpreter','latex','fontsize',14) ;
ylabel('$\tilde G_{xy}$ (MPa)','interpreter','latex') ;
xlabel('mode set','interpreter','latex') ;

sgtitle('Finnish Spruce Tonewood','interpreter','latex')





