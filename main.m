%---- test VK coefficient calc
clear all
close all
clc


%------------------------------------------------------------------------
% custom params
rho     = 8000 ;
E       = 2e11 ;
nu      = 0.3 ;
Lz      = 1e-3 ;
Lx      = 0.6 ;
Ly      = 0.4 ;
hfrac   = 0.01 ;   %-- computed as a fraction of sqrt(Lx*Ly)
Nmodes  = 5 ;
%BCs Transv
BCsPhi  = [0 0 ; 1e15 1e15 ; 0 0 ; 0 0] ;

%BCs Airy
BCsPsi  = [1e15 1e15 ; 1e15 1e15 ; 1e15 1e15 ; 1e15 1e15] ;
%-- NB: these represent mathematically "clamped" BCs, but "free" physically
%-- this choice enables the "triple self-adjointness" property, so best to keep this as is

Ntensor = Nmodes;
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% derived params
h       = sqrt(Lx*Ly)*hfrac ;
ldim    = [Lx Ly Lz] ;


[Om,Phi,Nx,Ny,~,~]       = magpie(rho,E,nu,ldim,h,BCsPhi,Nmodes,"none",1) ;
[~,Psi,~,~,~,zetafourth] = magpie(rho,E,nu,ldim,h,BCsPsi,Nmodes,"none",1) ;
zeta = (zetafourth).^(1/4) ;

H = zeros(Ntensor,Ntensor,Ntensor) ;
E = zeros(Ntensor,Ntensor,Ntensor) ;

Dxx = DxxBuild(Nx,Ny,h) ;
Dyy = DyyBuild(Nx,Ny,h) ;
Dx  = DxBuild(Nx,Ny,h) ;
Dy  = DyBuild(Nx,Ny,h) ;

for k = 1 : Ntensor
    Phik = Phi(:,k) ; Psik = Psi(:,k) ;
    for p = 1 : Ntensor
        Phip = Phi(:,p) ;
        for q = 1 : Ntensor
            Phiq = Phi(:,q) ; Psiq = Psi(:,q) ;

            LPhipPhiq = vkOperator(Phip,Phiq,Dx,Dy,Dxx,Dyy) ;
            LPhipPsiq = vkOperator(Phip,Psiq,Dx,Dy,Dxx,Dyy) ;

            H(k,p,q) = trapzIntcalc(Psik.*LPhipPhiq,h,Nx,Ny) ;
            E(k,p,q) = trapzIntcalc(Phik.*LPhipPsiq,h,Nx,Ny) ;

        end
    end
end
%%
for k=1:Nmodes
    for p=1:Nmodes
        for q=1:Nmodes
            symmetryH=H(k,p,q)-H(k,q,p)
            %symmetryE=H(k,p,q)-E(q,p,k)
        end
    end
end














