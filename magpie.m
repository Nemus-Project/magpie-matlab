%% MAGPIE Modal analysis of 2D Plates with Generalised Elastic Boundary Conditions
%   [Om,Q,Nx,Ny,biHarm,Dm] = MAGPIE (density, Youngs, poisson, dim, h, BCs, Number_of_modes, plot_type, normalisation)
%
%%   Arguments:
%       rho          %-- density [kg/m^3]
%       E            %-- Young's mod [Pa]
%       nu           %-- poisson's ratio
%
%       3 element array  representing  [x,y,z] dimensions of plate
%       ldim = [Lx,  %-- length along x [m]
%               Ly,  %-- length along y [m]
%               Lz]  %-- thickness [m]
%
%       h            %-- grid spacing
%
%       2 column array of eleastic boundary constants around each edge of
%       the plate.
%
%       Column 1 flexural constant
%       Column 2 rotational constant
%
%       BCs = [K0y, R0y;
%              Kx0, Rx0;
%              KLy, RLy;
%              KxL, RxL];
%
%       Nm      %-- number of modes to compute
%
%       plot_type %-- Select from 'chladni', '3d', 'none'
%
%       normalisation %-- boolean which dictates wether eigen vectors are
%       normalised
%
%% Returns:
%           Om      : Angular modal frequencies
%           Q       : A matrix of column eigenvector(s)
%           Nx      : Grid points along the x-axis
%           Ny      : Grid points along the y-axis
%           biHarm  : Biharmonic Matrix for the plate
%           Dm      : the eigenvalues of the biharmonic matrix
%
%% Example:
%
%           %% physical and elastic parameters
%           Lx = 0.10; Ly = 0.08; Lz = 0.81e-3;
%           ldim = [Lx Ly Lz];       % plate dimensions [x, y, z] in metres
%           E    = 1.01e+11 ;        %-- Young's mod [Pa]
%           rho  = 8765 ;            %-- density [kg/m^3]
%           nu   = 0.3 ;             %-- poisson's ratio
%           Nm   = 16;               %-- number of modes to compute
%           h    = sqrt(Lx*Ly)*0.01; %-- Grid Spacing
%           BCs = ones(4,2) * 1e15   %-- elastic constants around the edges
%
%           [Om,Q,Nx,Ny,biHarm,Dm] = magpie(rho, E, nu, ldim, h, BCs, Nm,'none');
%
function [Om,Q,Nx,Ny,biHarm,Dm] = magpie(rho,E,nu,ldim,h,BCs,Nm,plot_type,shouldNormalise)
%

%% Variable Arguments
if nargin < 9
    shouldNormalise = false;
end
if nargin < 8
    plot_type = 'none';
end
if nargin < 7
    Nm = 0;
end
%% Argument Validation
validateattributes(rho,  {'double'}, {'nonempty','positive'},'magpie.m','Density (rho)',1);
validateattributes(E,    {'double'}, {'nonempty','positive'},'magpie.m','Young`s Modulus (E)', 2);
validateattributes(nu,   {'double'}, {'nonempty','positive'},'magpie.m','Poisson Number (nu)',3);
validateattributes(ldim, {'double'}, {'numel', 3, 'positive'},'magpie.m','Plate Dimensions (ldim)',4);
validateattributes(h,    {'double'}, {'nonempty','positive'},'magpie.m','Grid Spacing (h)',5);
validateattributes(BCs,  {'double'}, {'size', [4,2], 'nonnegative'},'magpie.m','Boundary Conditions (BCs)',6);
validateattributes(Nm,   {'numeric'}, {'integer','nonnegative'},'magpie.m','Number of Modes (Nm)',7);
validatestring(plot_type,["chladni","3D","none"]);

%% Unpack array variables
pack_ldim = num2cell(ldim);
[Lx, Ly, Lz] = pack_ldim{:};

%%--- derived parameters (don't change here)
D = E * Lz^3 / 12 / (1-nu^2);
Nx      = floor(Lx/h) ;
Ny      = floor(Ly/h) ;

%% Build Biharmonic
biHarm = bhmat(BCs,[Nx Ny], h, Lz, E, nu);

%% EIGENVALUES

Nmodes = (Nx+1)*(Ny+1) ;
if Nm
    Nmodes = Nm ;
end
[Q,Dm] = eigs(biHarm,Nmodes,'smallestabs') ;
[~,indSort] = sort(diag((Dm))) ;
Q = Q(:,indSort) ;

Dm    = diag(Dm) ;
Om    = sqrt(abs(Dm))*sqrt(D/rho/Lz) ;

if shouldNormalise
    for nQ = 1 : Nmodes
        Qtemp   = Q(:,nQ) ;
        Qnorm   = trapzIntcalc(Qtemp.*Qtemp,h,Nx,Ny) ;
        Qtemp   = Qtemp / sqrt(Qnorm) ;
        Q(:,nQ) = Qtemp ;
    end
end

%% Plotting

switch plot_type
    case 'chladni'
        subs = ceil(sqrt(Nmodes));

        colormap('copper') ;
        cmp = colormap;
        cmp = flipud(cmp);
        colormap(cmp);

        for m = 1 : Nmodes
            mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
            subplot(subs,subs,m)
            mesh(3e3*real(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;
            view(2);
            axis equal;
            axis tight;
            axis off;
            clim([0.00005 0.002]);
        end

    case '3D'

        subs = ceil(sqrt(Nmodes));
        colormap('parula') ;
        xax = (0:Nx)*h ;
        yax = (0:Ny)*h ;
        [X,Y] = meshgrid(xax,yax) ;

        for m = 1 : Nmodes
            mdShape = reshape(Q(:,m),[(Ny+1),(Nx+1)]) ;
            subplot(subs,subs,m)
            mesh(X,Y,3000*(mdShape),(abs(mdShape)),'FaceColor','texturemap') ;
        end
end

end
