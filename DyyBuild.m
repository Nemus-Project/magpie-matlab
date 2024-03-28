function Dyy = DyyBuild(Nx,Ny,h)

%------------------------------------------------------------------------
% Dy
DyyBlk = sparse(Ny+1,Ny+1) ;


for n = 2 : Ny

    DyyBlk(n,n-1) = 1/h^2 ; DyyBlk(n,n+1) = 1/h^2  ; DyyBlk(n,n) = -2/h^2  ;

end

DyyBlk(1,1) = 2/h^2 ; DyyBlk(1,2) = -5/h^2 ; DyyBlk(1,3) = 4/h^2 ; DyyBlk(1,4) = -1/h^2 ;
DyyBlk(end,end) = 2/h^2 ; DyyBlk(end,end-1) = -5/h^2 ; DyyBlk(end,end-2) = 4/h^2 ; DyyBlk(end,end-3) = -1/h^2 ;

Dyy = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

for m = 1 : Nx + 1

    Dyy((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1)) = DyyBlk ;

end
%------------------------------------------------------------------------
% Dyy = (1/h^2)*Dyy
d0=-2*ones((Nx+1)*(Ny+1),1);
d0([1:Ny+1:end])=2; 
%Ny:Ny+1:end
end
