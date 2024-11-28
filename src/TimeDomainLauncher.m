%-------------------------------------------------------------------------
%             Time Domain Launcher 
%   Orthotropic Plates with Ribs, Static loads and Stiffeners 
%                Dr M Ducceschi
%             University of Bologna
%                  30 Sep 2024
%-------------------------------------------------------------------------

clear all
close all
clc

%-------------------------------------------------------------------------
% CUSTOM PARAMETERS (user can change these)

%--------------------
%-- general parameters
fmax = 5000 ; % maximum frequency to be computed
ppw  = 3 ; % points per wavelength at maximum frequency. Choose 3 <= ppw 
T    = 0.5 ; % total simulation time [s]
%--------------------

%--------------------
%-- plate parameters
rho   = 390 ;    % density [kg / m^3]
nux   = 0.39 ;   % poisson ratio
Ex    = 10.9e9 ; % young's mod along x [Pa]
Ey    = 0.64e9 ; % young's mod along y [Pa]
Gxy   = 0.58e9 ; % shear mod [Pa]
Lx    = 0.6 ;    % edge lenght x [m]
Ly    = 0.6 ;    % edge length y [m]
Lz    = 1.0e-3 ; % thickness [m]

%-- elastic constants around the edges
KRmat = [1e4,1e5; % Kx0 Rx0 
    1e4, 1e4;     % K0y R0y
    0e10, 0e10;   % KxL RxL
    0e10, 0e10] ; % KLy RLy
%--------------------

%--------------------
%-- input / output locations
x_p = [0.57,0.608] ; % frac of Lx Ly 
outMat = [0.57,0.86;
    0.57,0.27] ; % frac of Lx Ly, each row represents an output point
%--------------------

%--------------------
%-- forcing
forceType = 1 ; % 1 = impulse, 2 = sinusoid, 3 = input_file
if forceType == 1 
forceParams = [0.0007,5,0.5] ; % [time of contact (s), max Amplitude (N), noise modulation]
elseif forceType == 2
    forceParams = [100,50,0.2] ; % [frequency of sine (Hz), max Amplitude (N), noise modulation]
elseif forceType == 3
    forceParams = 'testInput.wav' ;
end
%--------------------

%--------------------
%-- loss parameters
dampVec = [0.3, 0.23, 500] ; %[t60_0Hz (s),t60_f1 (s), f1 (Hz)] NB you MUST use t60_f1 < t60_0Hz for stability
%--------------------


%--------------------
%-- braces parameters
Nribs = 4 ; % number of braces 

Eb = [10e9, 11e9, 10.2e9, 200e9].' ; % youngs moduli
Lzb = [3e-3,2e-3,2e-3,1e-3].' ; % thicknesses 
bb = [3e-2,2e-2,2e-2,1e-2].' ; % width cross section
rhob = [400,390,410,8000].' ; % densities

% rib coordinates along x (start and end) AS A FRACTION OF Lx
x_beam_coord = ...
    [0.2,0.8;
    0.2,0.8;
    0.2,0.8;
    0.2,0.8] ;

% rib coordinates along y (start and end) AS A FRACTION OF Ly
y_beam_coord = ...
    [0.1,0.2;
    0.2,0.3
    0.3,0.4
    0.4,0.5] ;
%--------------------

%--------------------
%-- static loads and stiffeners parameters
Nlump        = 2 ;                  % number of lumped elements
x_lump_coord = [0.2,0.4].'  ;       % x coordinates of lumped elements AS A FRACTION OF Lx
y_lump_coord = [0.7,0.87].' ;       % y coordinates of lumped elements AS A FRACTION OF Ly
Mlump        = [0.2,0.01].' ;  % masses [kg]
%--------------------

%--------------------
%-- plot parameters parameters
cmap = cmaps(2) ; % select colormap 1 = RedBlue, 2 = GreenPurple, 3 = OrangeGreen, 4 = OrangePurple
LivePlot = 1 ; % 1 : live plot on
RefreshRate = 5 ; % 1 = play all frames, 2 = play one out of two frames, etc
absPlot = 0 ; % 1 = plots absolute value (colormap will adjust accordingly)
FilmRec = 0 ; % 1= record video to file
%--------------------

% END CUSTOM PARAMETERS 
%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% DERIVED PARAMETERS (user cannot change these)

%-- plate grid spacings and grid points
Evec = [Ex,Ey,Gxy] ;
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

Lb = sqrt((x_beam_coord(:,1)-x_beam_coord(:,2)).^2 + (y_beam_coord(:,1)-y_beam_coord(:,2)).^2) ;
Nb = round(Lb./hb) ; % number of points
hb = Lb./Nb ; % actual grid spacings 

beamParams = [rhob,Eb,Ib,Ab,Nb,hb] ;
beamCoord = [x_beam_coord; y_beam_coord] ;


%-- lumped elements parameters
lumpCoord = [x_lump_coord; y_lump_coord] ;


%-- plot parameters
close all
xax = linspace(0,Lx,Nx+1) ;
yax = linspace(0,Ly,Ny+1) ;
[X,Y] = meshgrid(xax,yax) ;


%-------------------------------------------------------------------------


%-------------------------------------------------------------------------
% EIGENVALUE PROBLEM
[outs,stresses,fs,fin] = time_domain_sim(rho,Evec,nux,Lvec,hvec,Nvec,KRmat,fmax,T,x_p,outMat,Nribs,beamParams,beamCoord,Nlump,Mlump,lumpCoord,dampVec,forceType,forceParams,LivePlot,RefreshRate,absPlot,FilmRec,cmap) ;
%-------------------------------------------------------------------------



