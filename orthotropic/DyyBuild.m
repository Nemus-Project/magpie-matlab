function Dyy = DyyBuild(Nx,Ny,h,method)

%------------------------------------------------------------------------

if nargin < 4
    method = 'diag';
end

%% Validate Arguments
valid_methods = ["blk", "diag"];

validateattributes(Nx, {'numeric'}, {'integer','positive'});
validateattributes(Ny, {'numeric'}, {'integer','positive'});
validateattributes(h,  {'numeric'}, {'real','positive'});
method = validatestring (method,    valid_methods);
% Dy
switch(method)
    case('blk')
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
    case('diag')
        d0=-2*ones((Nx+1)*(Ny+1),1);
        d0([1:Ny+1:end,Ny+1:Ny+1:end])=2;

        d1 = ones(Ny+1,1);  %%Matthew's trick
        d1([1,end]) = [-5,0];
        d1 = repmat(d1,Nx+1,1);
        d1=d1(1:end-1);

        d2=zeros((Nx+1)*(Ny+1)-2,1);
        d2(1:Ny+1:end)=4;

        d3=zeros((Nx+1)*(Ny+1)-3,1);
        d3(1:Ny+1:end)=-1;


        % assert(all(diag(Dyy,0)-d0 ==0),'diag 0 is pure man facked')
        % assert(all(diag(Dyy,1)-d1 ==0),'diag 1 is pure man facked')
        % assert(all(diag(Dyy,2)-d2 ==0),'diag 2 is pure man facked')
        % assert(all(diag(Dyy,3)-d3 ==0),'diag 3 is pure man facked')
        % assert(all(diag(Dyy,-1)-flipud(d1) ==0),'diag -1 is pure man facked')

        d=[[flipud(d3);[0;0;0]], [flipud(d2);[0;0]],...
            [flipud(d1);[0]],...
            d0,...
            [[0];d1],...
            [[0;0];d2], [[0;0;0];d3]];

        %d=[[d3Ny1;zeros(3*(Ny+1),1)]]

        dN= (-3:3);
        %dN= (Ny+1)*(3);
        %
        Dyy = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
        Dyy = (1/h^2)*spdiags(d,dN,Dyy);

        %assert(all(all(Dyy-Dyy ==0)),'Dyytest is pure man facked')
end
end
