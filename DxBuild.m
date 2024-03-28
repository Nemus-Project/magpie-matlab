function Dx = DxBuild(Nx,Ny,h)

%------------------------------------------------------------------------
% Dx
Dx = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 2 : Nx

    Dx((m-1)*(Ny+1)+1:m*(Ny+1),(m-2)*(Ny+1)+1:(m-1)*(Ny+1)) = speye(Ny+1)*(-1/2/h) ;
    Dx((m-1)*(Ny+1)+1:m*(Ny+1),m*(Ny+1)+1:(m+1)*(Ny+1)) = speye(Ny+1)*(1/2/h) ;

end

Blk00 = speye(Ny+1,Ny+1)*(-3/2/h) ;
Blk01 = speye(Ny+1,Ny+1)*(2/h) ;
Blk02 = speye(Ny+1,Ny+1)*(-1/2/h) ;

Dx(1:Ny+1,1:Ny+1) = Blk00 ;
Dx(1:Ny+1,Ny+2:2*Ny+2) = Blk01 ;
Dx(1:Ny+1,2*Ny+3:3*Ny+3) = Blk02 ;

Dx(end-Ny:end,end-Ny:end) = -Blk00 ;
Dx(end-Ny:end,end-2*Ny-1:end-Ny-1) = -Blk01 ;
Dx(end-Ny:end,end-3*Ny-2:end-2*Ny-2) = -Blk02 ;
%------------------------------------------------------------------------

end
