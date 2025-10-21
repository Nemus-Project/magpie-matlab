%-------------------------------------------------------------------------
%             Frequency Domain Launcher 
%   Orthotropic Plates with Ribs, Static loads and Stiffeners 
%                Dr M Ducceschi
%             University of Bologna
%                  30 Sep 2024
%-------------------------------------------------------------------------

clear all
close all
% clc

%-------------------------------------------------------------------------
% CUSTOM PARAMETERS (user can change these)

%--------------------
%-- general parameters
fmax = 5000 ; % maximum frequency to be computed
ppw  = 50 ; % points per wavelength at maximum frequency. Choose 3 <= ppw 
Nmodes = 9 ; % select total number of modes to be computed. If Nmodes = [], all available modes are computed
%--------------------

%--------------------
%-- plate parameters
rho   = 457 ;    % density [kg / m^3]
Lx    = 0.212 ;  % edge lenght x [m]
Ly    = 0.108 ;  % edge length y [m]
Lz    = 0.0045 ; % thickness [m]

%-- guessed values
Ex0    = 13.1e9 ;  % young's mod along x [Pa]
Ey0    = 0.881e9 ; % young's mod along y [Pa]
Gxy0   = 0.504e9 ; % shear mod [Pa]
nux0   = 0.4 ;     % poisson ratio


%-- elastic constants around the edges
KRmat = [0e10,0e10; % Kx0 Rx0 
    1e10, 1e10;     % K0y R0y
    0e10, 0e10;   % KxL RxL
    0e10, 0e10] ; % KLy RLy
%--------------------

%--------------------
%-- braces parameters
Nribs = 0 ; % number of braces 

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
Nlump        = 0 ;                  % number of lumped elements
x_lump_coord = [0.2,0.4].'  ;       % x coordinates of lumped elements AS A FRACTION OF Lx
y_lump_coord = [0.7,0.87].' ;       % y coordinates of lumped elements AS A FRACTION OF Ly
Mlump        = [0.2,0.01].' ;  % masses [kg]
%--------------------

%--------------------
%-- plot parameters parameters
cmap = cmaps(4) ; % select colormap 1 = RedBlue, 2 = GreenPurple, 3 = OrangeGreen, 4 = PurpleOrange
NN = 1 ; % first mode number to be plotted 
Nplots = 9 ; % select 3,6 or 9. If another number is selected, it is defaulted to 3. Displayed plots are NN + (0:Nplots - 1)
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
Nvec = [Nx,Ny];  % grid points
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
[Om,Q,K,M,Qtot,OmTot] = freq_domain_sim(rho,Evec,nux0,Lvec,hvec,Nvec,KRmat,fmax,Nmodes,Nribs,beamParams,beamCoord,Nlump,Mlump,lumpCoord) ;
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% PLOT RESULTS
modal_plotter(cmap,Nplots,NN,X,Y,Q,Lvec,Nvec,Nribs,Nlump,beamParams,beamCoord,lumpCoord);
%-------------------------------------------------------------------------

ModeFrequencies = Om/2/pi;
digits(3);
ModeFrequencies = vpa(ModeFrequencies)
