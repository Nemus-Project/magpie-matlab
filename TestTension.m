clear all
%close all
clc
%%
% This code computes the modes of the displacement and Airy's stress of a
% nonlinear Föpple-von Karmàn plate, and the coupling coefficients between
% the modes. It takes as an input:
%
%   PHYSICAL PARAMETERS
%       E - Young modulus                           | float
%       rho - Density                               | float
%       nu - Poisson's ratio                        | float
%       Lx - x length of the plate                  | float
%       Ly - y length of the plate                  | float
%       Lz - Thickness of the plate                 | float
%       BCsPhi - Displacement boundary conditions   | 4 x 2 vector
%       BCsPsi - Stress boundary conditions         | 4 x 2 vector
%
%   NUMERICAL PARAMETERS
%       
%       Nx - Number of points in the x direction    | float
%       Nvec - Vector containing a test Nx and Nx   | 1 x 2 vector
%       Nmodes  - Number of modes                   | float
%
% The code then uses the magpie submodule to compute the linear modes of
% the plate for both displacement and stress, compare them with a
% refference through the functions, eigenmac and eigensign, and compute the
% nonlinear coupling coefficient through vkoperator.
%
% It outputs :
%
%       Hv - Coupling coeficients                   | Nmodes x Nmodes x Nmodes matrix
%       Om - Angular frequency of the Phi modes     | Nmodes x 1 vector
%       Om2 - Angular frequency of the Psi modes    | Nmodes x 1 vector
%       zetafourth - Eigenvalues of the stress      | Nmodes x 1 vector
%       Phi - Eigenvectors of displacement          | (Nx+1*Ny+1) x Nmodes matrix
%       Psi - Eigenvectors of stress                | (Nx+1*Ny+1) x Nmodes matrix
%
% The last section saves all the important variables in a file located in
% the 'param' folder, for later use by the 'main.m' code.


%% Adding path to the submodule magpie

addpath ./private/magpie

%% Variable declaration

% Physical parameters
E       = 4e+9 ;
rho     = 900 ; 
nu      = 0.4 ;
%sigma   = rho*2e-5;
Tvec       = [0:16]*0.2250;
freqshift=zeros(17,4)
Lzvec=[2e-5,4e-5,8e-5];
Lxvec=[4.3e-2,4e-2,3.8e-2,3.5e-2]
for it=1:4
Lz      = 2e-5 ;
Lx      = Lxvec(it);
Ly      = Lx;
D       = E * Lz^3 / 12 / (1-nu^2);
%sigma   = rho*Lz;

for itit = 1:17
T       = Tvec(itit)


% E       =  200000000000.0 ;
% rho     = 8000 ; 
% nu      = 0.3 ;
% Lz      = 0.001 ;
% Lx      =  0.2 ;
% Ly      = 0.4 ;
% D       = E * Lz^3 / 12 / (1-nu^2);
% T       = 0*D;

% Numerical parameters
Nmodes  = 5 ;
Nx=400;

% Derived values
Nvec=[100 Nx];
npts=length(Nvec);


% BCs Displacement
BCsPhi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;

%BCs Airy's stress
BCsPsi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%-- NB: these represent mathematically "clamped" BCs, but "free" physically
%-- this choice enables the "triple self-adjointness" property, so best to keep this as is

ldim    = [Lx Ly Lz] ;
%%
for iter=1:npts

    h=Lx/Nvec(iter); % Defining h

    [Om,Phi,Nx,Ny,~,~] = magpie(rho,E,nu,T,ldim,h,BCsPhi,Nmodes,"none",true) ; % see magpie doc

    if iter ==1

        Phiref=Phi; % Setting a refference for Psi
        Nxref=Nx; % Setting a refference for Nx
        Nyref=Ny; % Setting a refference for Ny

    else

        [Phi,Om] = eigenMAC(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly,Om); % Ensures that mode order is consistent

        Phi = eigensign(Phiref,Nxref,Nyref,Phi,Nx,Ny,h,Nmodes,Lx,Ly); % Ensures that the polarization is consistent

    end

    [Om2,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,T,ldim,h,BCsPsi,Nmodes,"none",true) ;% see magpie doc

    if iter ==1

        Psiref=Psi; % Setting a refference for Psi

    else

        [Psi,Om2] = eigenMAC(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly,Om2);% Ensures that mode order is consistent

        Psi = eigensign(Psiref,Nxref,Nyref,Psi,Nx,Ny,h,Nmodes,Lx,Ly);% Ensures that the polarization is consistent

    end

    disp(iter)

    zeta = (zetafourth).^(1/4) ;

    Hv = zeros(Nmodes,Nmodes,Nmodes) ;

    Ev = zeros(Nmodes,Nmodes,Nmodes) ;

    Dxx = DxxBuild(Nx,Ny,h) ;
    Dyy = DyyBuild(Nx,Ny,h) ;
    Dxy = DxyBuild(Nx,Ny,h) ;

    tic;
    for k = 1 : Nmodes

        Phik = Phi(:,k) ; Psik = Psi(:,k) ;

        for p = 1 : Nmodes
            
            Phip = Phi(:,p);

            for q = p : Nmodes

                Phiq = Phi(:,q) ; Psiq = Psi(:,q);

                LPhipPhiq = vkOperator(Phip,Phiq,Dxy,Dxx,Dyy) ;

                %norm_k = trapzIntcalc(Psik.*Psik,h,Nx,Ny);
                
                Hv(k,q,p) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny);

                Hv(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny); %Coupling coefficient tensor
            end
        end
    end
end
%% Defining the numerical space
xax = (0:Nx)*h ;
yax = (0:Ny)*h ;
[X,Y] = meshgrid(xax,yax) ;

% %% Save parameters
% 
% if ~exist("./param/", 'dir')
%     mkdir("./param/")
% end
% 
% filename='Purple20modes';
% 
% save(['./param/' filename '.mat'],'rho','E','nu','Lz','Lx','Ly','Nmodes','Phi','Om','Psi','Om2','Nx','Ny','h','X','Y','zetafourth','BCsPhi','BCsPsi','Hv');

freqshift(itit,it)=Om(1)./(2*pi)
end
end
%%
% newcolors = [
%              0.8 0.4 0
%              0.69 0.28 0.82
%              0 0.7 0
%              0 0.4470 0.7410
%              0 0 0
%              0.1 0.50 0.95
%              0 0 0 ];
newcolors = [
             0.8 0.4 0
             0 0 0
             0.8 0.4 0
             0 0 0
             ];
figure
colororder(newcolors)
plot(freqshift(:,1:2),Tvec,"LineWidth",3,"Marker","o","MarkerSize",15)
hold on
plot(freqshift(:,3:4),Tvec,"LineWidth",3,"Linestyle","--","Marker","o","MarkerSize",15)
set(gca,"FontSize",26)
ylabel("T_0 (N/m)")
xlabel("Frequency of the first mode (Hz)")