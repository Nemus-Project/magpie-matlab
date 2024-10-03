function Dxy = DxyBuild(Nx,Ny,hx,hy,method)
% Dxx = DXXBUILD(Nx,Ny,h,method)
%
%       Nx     : number of grid points x-axis
%       Ny     : number of grid points y-axis
%       hx     : grid spacing along x
%       hy     : grid spacing along y
%       method : method for building matrix: 'blk' or 'diag'
%
%------------------------------------------------------------------------

% Parse args
if nargin < 4
    method = 'diag';
end
%% Validate Arguments
valid_methods = ["blk", "diag"];

validateattributes(Nx, {'numeric'}, {'integer','positive'});
validateattributes(Ny, {'numeric'}, {'integer','positive'});
%validateattributes(h,  {'numeric'}, {'real','positive'});
method = validatestring (method,    valid_methods);
% Validate args

% Dy
switch(method)
    case('blk')


        DyBlk = sparse(Ny+1,Ny+1) ;


        for n = 2 : Ny

            DyBlk(n,n-1) = -1/2/hy ; DyBlk(n,n+1) = 1/2/hy  ;

        end

        DyBlk(1,1) = -3/2/hy ; DyBlk(1,2) = 2/hy ; DyBlk(1,3) = -1/2/hy ;
        DyBlk(end,end) = 3/2/hy ; DyBlk(end,end-1) = -2/hy ; DyBlk(end,end-2) = 1/2/hy ;

        Dy = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

        for m = 1 : Nx + 1

            Dy((m-1)*(Ny+1)+1:m*(Ny+1),(m-1)*(Ny+1)+1:m*(Ny+1)) = DyBlk ;

        end

        % Dx
        Dx = sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;

        for m = 2 : Nx

            Dx((m-1)*(Ny+1)+1:m*(Ny+1),(m-2)*(Ny+1)+1:(m-1)*(Ny+1)) = speye(Ny+1)*(-1/2/hx) ;
            Dx((m-1)*(Ny+1)+1:m*(Ny+1),m*(Ny+1)+1:(m+1)*(Ny+1)) = speye(Ny+1)*(1/2/hx) ;

        end

        Blk00 = speye(Ny+1,Ny+1)*(-3/2/hx) ;
        Blk01 = speye(Ny+1,Ny+1)*(2/hx) ;
        Blk02 = speye(Ny+1,Ny+1)*(-1/2/hx) ;

        Dx(1:Ny+1,1:Ny+1) = Blk00 ;
        Dx(1:Ny+1,Ny+2:2*Ny+2) = Blk01 ;
        Dx(1:Ny+1,2*Ny+3:3*Ny+3) = Blk02 ;

        Dx(end-Ny:end,end-Ny:end) = -Blk00 ;
        Dx(end-Ny:end,end-2*Ny-1:end-Ny-1) = -Blk01 ;
        Dx(end-Ny:end,end-3*Ny-2:end-2*Ny-2) = -Blk02 ;

        Dxy=Dx*Dy;



        %-------
    case('diag')
        d0=zeros((Nx+1)*(Ny+1),1);
        d0([1,end])=2.2500;
        d0([Ny+1,end-Ny])=-2.2500; %%central diagonal

        d1=zeros((Nx+1)*(Ny+1)-1,1);
        d1(2:Ny)=-0.7500;
        d1(end-Ny+2:end)=0.7500;
        d1(1)=-3.0000;
        d1(end-Ny+1)=3.0000; %% 1 off diagonal

        d2=zeros((Nx+1)*(Ny+1)-2,1);
        d2(1)=0.7500;
        d2(end-Ny+2)=-0.7500;  %% 2 off diagonal

        dNym1=zeros((Nx+1)*(Ny+1)-(Ny-1),1);
        dNym1(Ny+1)=1.0000;
        dNym1(2*Ny+2:Ny+1:end)=0.2500; % Ny-1 diagonal

        dNy=zeros((Nx+1)*(Ny+1)- Ny,1);
        mot = [0;-0.2500*ones(Ny-1,1);-1];
        dNy = [repmat(mot,Nx,1);0];
        dNy(2:Ny)=-1.0000;
        dNy(Ny+1)=-4.0000;% Ny diagonal

        dNyp1=zeros((Nx+1)*(Ny+1)-(Ny+1),1);
        dNyp1(1)=-3.0000;
        dNyp1(Ny+1)=3.0000;
        dNyp1(2*Ny+2:Ny+1:end)=0.7500;
        dNyp1(Ny+2:Ny+1:end)=-0.7500; % Ny+1 diagonal

        dNyp2=zeros((Nx+1)*(Ny+1)- (Ny+2),1);
        mot = [1;0.2500*ones(Ny-1,1);0];
        dNyp2 = repmat(mot,Nx,1);
        dNyp2 = dNyp2(1:end-1);
        dNyp2(2:Ny)=1.0000;
        dNyp2(1)=4.0000;% Ny+2 diagonal


        dNyp3=zeros((Nx+1)*(Ny+1)-(Ny+3),1);
        dNyp3(1)=-1.0000;
        dNyp3(Ny+2:Ny+1:end)=-0.2500; % Ny+3 diagonal

        d2Ny=zeros((Nx+1)*(Ny+1)-(2*Ny),1);
        d2Ny(Ny+1)=-0.2500;% 2*Ny diagonal
        

        d2Nyp1=zeros((Nx+1)*(Ny+1)-(2*Ny+1),1);
        d2Nyp1(2:Ny)=0.2500;
        d2Nyp1(Ny+1)=1.0000;% 2*Ny+1 diagonal

        d2Nyp2=zeros((Nx+1)*(Ny+1)-(2*Ny+2),1);
        d2Nyp2(1)=0.7500;
        d2Nyp2(Ny+1)=-0.7500;% 2*Ny+2 diagonal

        d2Nyp3=zeros((Nx+1)*(Ny+1)-(2*Ny+3),1);
        d2Nyp3(1)=-1.0000;
        d2Nyp3(2:Ny)=-0.2500;% 2*Ny+3 diagonal

        d2Nyp4=zeros((Nx+1)*(Ny+1)-(2*Ny+4),1);
        d2Nyp4(1)=0.2500;% 2*Ny+4 diagonal







        %  assert(all(diag(Dxy,0)-d0 ==0),'diag 0 is pure man facked')
        % 
        % assert(all(diag(Dxy,1)-d1 ==0),'diag 1 is pure man facked')
        % 
        % assert(all(diag(Dxy,2)-d2 ==0),'diag 2 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny-1)-dNym1 ==0),'diag Ny-1 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny)-dNy ==0),'diag Ny is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny+1)-dNyp1 ==0),'diag Ny+1 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny+2)-dNyp2 ==0),'diag Ny+2 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny+3)-dNyp3 ==0),'diag Ny+3 is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny+1)-d2Nyp1 ==0),'diag 2*Ny+1 is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny+2)-d2Nyp2 ==0),'diag 2*Ny+2 is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny+3)-d2Nyp3 ==0),'diag 2*Ny+3 is pure man facked')
        % 
         % assert(all(diag(Dxy,2*Ny+4)-d2Nyp4 ==0),'diag 2*Ny+4 is pure man facked')


       

        d=[[flipud(d2Nyp4);zeros(1*(2*Ny+4),1)],...
            [flipud(d2Nyp3);zeros(1*(2*Ny+3),1)],...
            [flipud(d2Nyp2);zeros(1*(2*Ny+2),1)],...
            [flipud(d2Nyp1);zeros(1*(2*Ny+1),1)],...
            [flipud(d2Ny);zeros(1*(2*Ny),1)],...
            [flipud(dNyp3);zeros(1*(Ny+3),1)],...
            [flipud(dNyp2);zeros(1*(Ny+2),1)],...
            [flipud(dNyp1);zeros(1*(Ny+1),1)],...
            [flipud(dNy);zeros(1*(Ny),1)],...
            [flipud(dNym1);zeros(1*(Ny-1),1)],...
            [flipud(d2);0;0],...
            [flipud(d1);0],...
            d0,...
            [0;d1],...
            [0;0;d2],...
            [zeros(1*(Ny-1),1);dNym1],...
            [zeros(1*(Ny),1);dNy],...
            [zeros(1*(Ny+1),1);dNyp1],...
            [zeros(1*(Ny+2),1);dNyp2],...
            [zeros(1*(Ny+3),1);dNyp3],...
            [zeros(1*(2*Ny),1);d2Ny],...
            [zeros(1*(2*Ny+1),1);d2Nyp1],...
            [zeros(1*(2*Ny+2),1);d2Nyp2],...
            [zeros(1*(2*Ny+3),1);d2Nyp3],...
            [zeros(1*(2*Ny+4),1);d2Nyp4]];




        dN= [-(2*Ny+4),-(2*Ny+3),-(2*Ny+2),-(2*Ny+1),-(2*Ny),...
            -(Ny+3),-(Ny+2),-(Ny+1),-Ny,-(Ny-1),...
            -2,-1,0,1,2,...
            Ny-1, Ny, Ny+1,Ny+2,Ny+3,...
            2*Ny,2*Ny+1,2*Ny+2,2*Ny+3,2*Ny+4];


        %
        Dxy= sparse((Nx+1)*(Ny+1),(Nx+1)*(Ny+1)) ;
        Dxy = (1/h^2)*spdiags(d,dN,Dxy);

        % assert(all(diag(Dxy,0)-diag(Dxytest,0) ==0),'diag 0 is pure man facked')
        % 
        % assert(all(diag(Dxy,1)-diag(Dxytest,1) ==0),'diag 1 is pure man facked')
        % 
        % assert(all(diag(Dxy,2)-diag(Dxytest,2) ==0),'diag 2 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny-1)-diag(Dxytest,Ny-1) ==0),'diag Ny-1 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny)-diag(Dxytest,Ny) ==0),'diag Ny is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny+1)-diag(Dxytest,Ny+1) ==0),'diag Ny+1 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny+2)-diag(Dxytest,Ny+2) ==0),'diag Ny+2 is pure man facked')
        % 
        % assert(all(diag(Dxy,Ny+3)-diag(Dxytest,Ny+3) ==0),'diag Ny+3 is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny)-diag(Dxytest,2*Ny) ==0),'diag 2*Ny is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny+1)-diag(Dxytest,2*Ny+1) ==0),'diag 2*Ny+1 is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny+2)-diag(Dxytest,2*Ny+2) ==0),'diag 2*Ny+2 is pure man facked')
        % 
        % assert(all(diag(Dxy,2*Ny+3)-diag(Dxytest,2*Ny+3) ==0),'diag 2*Ny+3 is pure man facked')
        % 
        %  assert(all(diag(Dxy,2*Ny+4)-diag(Dxytest,2*Ny+4) ==0),'diag 2*Ny+4 is pure man facked')
        % 
        % assert(all(all(Dxy-Dxytest ==0)),'Dxytest is pure man facked')
end
end
