clear all
close all
clc

%% Biharm case stuff

%-- elastic constants around the edges
K0y     = 1e15 ;
R0y     = 0e15 ;
Kx0     = 1e15 ;
Rx0     = 0e15 ;
KLy     = 1e15 ;
RLy     = 0e15;
KxL     = 1e15 ;
RxL     = 0e15 ;

BCs=zeros(4,2);
BCs(:,1)=1e15;



%--- derived parameters (don't change here)
Lz=1e-3;
nu=0.3;
E=1e9;
D       = E * Lz^3 / 12 / (1-nu^2) ;
%k       =1/44100;
%kappa   =sqrt(D/(rho*Lz));
%h       =2*sqrt(kappa*k)
%h       = 0.67*(sqrt(Lx*Ly)*0.01);
%Nx      = floor(Lx/h) ;
%Ny      = floor(Ly/h) ;

%bhmat(BCs,[6,7],0.01,Lz,E,nu,'diag');

%%

nbpts=150;
blkt=zeros(nbpts,1);
olddiagt=zeros(nbpts,1);
diagt=zeros(nbpts,1);
nxny=zeros(nbpts,1);
for n=nbpts:-1:1
    Nx=10*n;
    Ny=1000;
    h=0.01;
    nxny(n)= Nx*Ny;
    tic
    %Dxy=DxyBuild(Nx,Ny,h,'blk');
    bhmat(BCs,[Nx,Ny],h,Lz,E,nu,'blk');
    blkt(n)=toc;

     tic
     
    %Dxy=DxyBuild(Nx,Ny,h,'diag');
    bhmat(BCs,[Nx,Ny],h,Lz,E,nu,'diag');
     
    olddiagt(n)=toc;
    

    % tic
    % for ite = 1:20
    % %Dyy=DyyBuild(Nx,Ny,h,'diag');
    % bhmat(BCs,[Nx,Ny],h,Lz,E,nu,'diag');
    % end
    % diagt(n)=toc;
    disp(n)


end
%%
figure
plot(nxny,blkt,nxny,olddiagt,Linewidth=3)
legend("Block","Diag")
xlabel("Number of operator points")
ylabel("Time (s)")
set(gca,'Fontsize',20)