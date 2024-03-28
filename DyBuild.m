function Dy = DyBuild(Nx,Ny,h)

%------------------------------------------------------------------------
% Dy
DyBlk = sparse(Ny+1,Ny+1) ;


for n = 2 : Ny

    DyBlk(n,n-1) = -1/2/h ; DyBlk(n,n+1) = 1/2/h  ;

end

DyBlk(1,1) = -3/2/h ; DyBlk(1,2) = 2/h ; DyBlk(1,3) = -1/2/h ;
DyBlk(end,end) = 3/2/h ; DyBlk(end,end-1) = -2/h ; DyBlk(end,end-2) = 1/2/h ;

Dy = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 1 : Nx + 1

    Dy((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1)) = DyBlk ;

end
%------------------------------------------------------------------------


end
