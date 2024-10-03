function J = spreadinterp(x_p,Nx,Ny,hx,hy,type)

%----------------------------------------------------------------------------------------
%  x_p  =  input point. Must be a two-element vector. Each element gives the fractional
%          location along x and y. Example x_p = [0.7,0.4] is equal to 0.7*Lx, 0.4*Ly in physical coordinates.
%          clearly, this must be satisifed 0 \leq x_p(1 or 2) \leq 1
%
%  Nx   = number of grid subintervals in the x direction

%  Ny   = number of grid subintervals in the y direction

%  hx   = grid spacing in the x-direction

%  hy   = grid spacing in the y-direction

%  type = 1: spreading (i.e., approximate a Dirac delta), 2: interpolator


nx     = floor(Nx*x_p(1)) ;
alx    = Nx*x_p(1) - nx ;
ny     = floor(Ny*x_p(2)) ;
aly    = Ny*x_p(2) - ny ;

J = zeros(Ny+1,Nx+1) ;

if nx < Nx && ny < Ny
    J(ny+1,nx+1)   = (1-alx)*(1-aly) ;
    J(ny+2,nx+1)   = alx*(1-aly) ;
    J(ny+1,nx+2)   = (1-alx)*aly ;
    J(ny+2,nx+2)   = alx*aly ;
elseif nx == Nx && ny < Ny
    J(ny+1,Nx+1)     = 1-aly ;
    J(ny+2,Nx+1)   = aly ;
elseif nx < Nx && ny == Ny
    J(Ny+1,nx+1)     = 1-alx ;
    J(Ny+1,nx+2)     = alx ;
else
    J(Ny+1,Nx+1) = 1 ;
end

J = reshape(J,(Nx+1)*(Ny+1),1) ;

if type == 1
    J = J/hx/hy ;
elseif type == 2
    J = J' ;
end

end

