%% TRAPZINTCALC
% short description
%% Syntax
%
% |trapzInt = trapzIntcalc(F,h,Nx,Ny)|
%
%% Description
%
% |trapzInt = trapzIntcalc(F,h,Nx,Ny)| returns....
%
%% Example
%
%
%
%% Input Arguments
%
% * |F|   : ...
% * |h|   : Grid spacing
% * |Nxy| : 2-element array of grid points |[Nx Ny]| used for calculation. 
%
function trapzInt = trapzIntcalc(F,h,Nxy)

Nx = Nx(1);
Ny = Ny(2);
%-----------------------------------------------------------
%--- define vector of weights
wVec = zeros(1,(Ny+1)*(Nx+1)) ;

wInterior       = ones(1,Ny+1) ;
wInterior(1)    = 0.5 ;
wInterior(end)  = 0.5 ;
wEdge           = 0.5 * ones(1,Ny+1) ;
wEdge(1)        = 0.25  ;
wEdge(end)      = 0.25  ;

for m = 2 : Nx
    wVec((m-1)*(Ny+1)+1:m*(Ny+1)) = wInterior ;
end

wVec(1,1:Ny+1)       = wEdge ;
wVec(1,end-Ny:end)   = wEdge ;

wVec = wVec * h^2 ;
%-----------------------------------------------------------

%-----------------------------------------------------------
%--- trapezoidal integration
trapzInt = wVec * F ;

end