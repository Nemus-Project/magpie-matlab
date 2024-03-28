function Dxx = DxxBuild(Nx,Ny,h,method)

%------------------------------------------------------------------------
% Dx
switch(method)
    case('blk')


        Dxx = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

        for m = 2 : Nx

            Dxx((m-1)*(Ny+1)+1:m*(Ny+1),(m-2)*(Ny+1)+1:(m-1)*(Ny+1)) = speye(Ny+1)*(1/h^2) ;
            Dxx((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1))     = speye(Ny+1)*(-2/h^2) ;
            Dxx((m-1)*(Ny+1)+1:m*(Ny+1),m*(Ny+1)+1:(m+1)*(Ny+1))     = speye(Ny+1)*(1/h^2) ;

        end

        Blk00 = speye(Ny+1,Ny+1)*(2/h^2) ;
        Blk01 = speye(Ny+1,Ny+1)*(-5/h^2) ;
        Blk02 = speye(Ny+1,Ny+1)*(4/h^2) ;
        Blk03 = speye(Ny+1,Ny+1)*(-1/h^2) ;

        Dxx(1:Ny+1,1:Ny+1) = Blk00 ;
        Dxx(1:Ny+1,Ny+2:2*Ny+2) = Blk01 ;
        Dxx(1:Ny+1,2*Ny+3:3*Ny+3) = Blk02 ;
        Dxx(1:Ny+1,3*Ny+4:4*Ny+4) = Blk03 ;

        Dxx(end-Ny:end,end-Ny:end) = Blk00 ;
        Dxx(end-Ny:end,end-2*Ny-1:end-Ny-1) = Blk01 ;
        Dxx(end-Ny:end,end-3*Ny-2:end-2*Ny-2) = Blk02 ;
        Dxx(end-Ny:end,end-4*Ny-3:end-3*Ny-3) = Blk03 ;
        %------------------------------------------------------------------------
        % Dxx = (1/h^2)*Dxx
    case('diag')
        d0=-2*ones((Nx+1)*(Ny+1),1);
        d0([1:Ny+1,(end-Ny):end])=2; %%TOTALLY ACCURATE VERSION

        dNy1=ones(Nx*(Ny+1),1);
        dNy1(1:(Ny+1))=-5;

        d2Ny1=zeros((Ny+1)*(Nx-1),1);
        d2Ny1(1:Ny+1)=4;

        d3Ny1=zeros((Ny+1)*(Nx-2),1);
        d3Ny1(1:Ny+1)=-1;

        %assert(all(diag(Dxx,0)-d0 ==0),'diag 0 is pure man facked')

        %assert(all(diag(Dxx,(Ny+1))-dNy1 ==0),'diag Ny+1 is pure man facked')

        %assert(all(diag(Dxx,2*(Ny+1))-d2Ny1 ==0),'diag 2*(Ny+1) is pure man facked')

        %assert(all(diag(Dxx,3*(Ny+1))-d3Ny1 ==0),'diag 3*(Ny+1) is pure man facked')


        d=[[flipud(d3Ny1);zeros(3*(Ny+1),1)], [flipud(d2Ny1);zeros(2*(Ny+1),1)],...
            [flipud(dNy1);zeros(1*(Ny+1),1)],...
            d0,...
            [zeros(1*(Ny+1),1);dNy1],...
            [zeros(2*(Ny+1),1);d2Ny1], [zeros(3*(Ny+1),1);d3Ny1]];

        %d=[[d3Ny1;zeros(3*(Ny+1),1)]]

        dN= (Ny+1)*(-3:3);
        %dN= (Ny+1)*(3);
        %
        Dxx = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
        Dxx = (1/h^2)*spdiags(d,dN,Dxx);

        %assert(all(all(Dxx-Dxx ==0)),'Dxxtest is pure man facked')
end
end
