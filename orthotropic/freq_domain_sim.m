function [Om,Q] = freq_domain_sim(rho,Evec,nux,Lvec,hvec,Nvec,KRmat,fmax,Nmodes,Nribs,beamParams,beamCoord,Nlump,lumpParams,lumpCoord)

%-----------------------------------------------------------------------
%-- unpack constants
Lz      = Lvec(3) ;
hx      = hvec(1) ;
hy      = hvec(2) ;
Nx      = Nvec(1) ;
Ny      = Nvec(2) ;
%-----------------------------------------------------------------------

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
    Klump = lumpParams(:,1) ;
    Mlump = lumpParams(:,2) ;
    for nL = 1 : Nlump

        x_lump = [lumpCoord(nL) ; lumpCoord(nL+Nlump)] ;
        LumpDist = [LumpDist,spreadinterp(x_lump,Nx,Ny,hx,hy,1)] ;

    end

    LumpSpread = hx*hy*LumpDist.' ;


    K = [K+LumpDist*diag(Klump)*LumpSpread, -LumpDist*diag(Klump); -diag(Klump)*LumpSpread, diag(Klump)] ;
    K = sparse(K) ;
    M = [M,zeros((Ny+1)*(Nx+1),Nlump); zeros(Nlump,(Ny+1)*(Nx+1)), sparse(diag(Mlump))] ;  M = sparse(M) ;
end


if isempty(Nmodes)
    Nmodes = (Nx+1)*(Ny+1) ;
end

[Q,Omsq] = eigs(K,M,Nmodes,'smallestabs') ;
[~,indSort] = sort(diag((Omsq))) ;
Q = Q(:,indSort) ;
Om = sqrt(abs(diag(Omsq))) ;



% keep values up to fmax
fCur = 0 ;
indCur = 0 ;
while fCur < fmax && indCur < Nmodes
    indCur = indCur + 1 ;
    fCur = Om(indCur) / 2 / pi ;
end
Om = Om(1:indCur) ;
Om(end)/2/pi
Q = Q(:,1:indCur) ;


%-----------------------------------------------------------------------


end



