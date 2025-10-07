%% BHCOEFS Biharmonic Coefficients for a given diagonal band
%   varargout = BHCOEFS(BCs,h,D,nu,order)
%   get individual coefficients of the isotropic biharmonic with 
%   elastic boundary constants
%
%% Arguments
%       2 column array of eleastic boundary constants around each edge of
%       the plate.
%
%       Column 1 flexural constant
%       Column 2 rotational constant
%
%      - BCs = [K0y, R0y;
%               Kx0, Rx0;
%               KLy, RLy;
%               KxL, RxL];
%      - h           %-- Grid spacing
%      - D           %-- Stiffness constant
%      - nu          %-- Poisson's ratio
%      - order       %-- Index of coefficient
%
%       valid_orders = ["D00", "D01", "D02", ...
%                       "D10", "D11", "D12", ...
%                       "D20", "D21", "D22"];
%
%% Returns
%         -  d0  : diagonal coefficient for index 0
%         -  d1  : diagonal coefficient for index 1
%         -  dc  : diagonal coefficient for central points
%         -  dM1 : diagonal coefficient for Ny - 1
%         -  dM  : diagonal coefficient for Ny
%
%% Example  
%           %% dm1   
%           
%           [d0,d1,dc,dM1,dM] = biharmdiag(BCs, h, D, nu,'-1');
%           
%           dm10 = d0(2)*a1;
%           dm10([1,Ny-1,Ny]) = [d0(1),d0(end-1),d0(end)];
%           dm11 = d1(2)*a1;
%           dm11([1,Ny-1,Ny]) = [d1(1),d1(end-1),d1(end)];
%           dm12 = dc(2)*a1;
%           dm12([1,Ny-1,Ny]) = [dc(1),dc(end-1),dc(end)];
%           dm1M1 = dM1(2)*a1;
%           dm1M1([1,Ny-1,Ny]) = [dM1(1),dM1(end-1),dM1(end)];
%           dm1M = dM(2)*a1;
%           dm1M([1,Ny-1,Ny]) = [dM(1),dM(end-1),dM(end)];
%
function varargout = bhcoefs(BCs,h,D,nu,order)
    %% Validate args
    
    validateattributes    (nu,     {'double'}, {'nonempty'});
    validateattributes    (D,      {'double'}, {'nonempty'});    
    validateattributes    (h,      {'double'}, {'nonempty'});    
    validateattributes    (BCs,    {'double'}, {'size', [4,2]});
    valid_orders = ["D00", "D01", "D02", ...
                    "D10", "D11", "D12", ...
                    "D20", "D21", "D22"];
    order = validatestring(order,    valid_orders);

    %% Unpack array variables
    pack_BCs = num2cell(BCs);
    [K0y, Kx0, KLy, KxL, R0y, Rx0, RLy, RxL] = pack_BCs{:};
    
    
    switch order
        case "D00"
            u00_c = ((2*K0y*Rx0*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*Rx0*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 12*Rx0*D^2*R0y*h^2*nu^2 + 24*Rx0*D^2*R0y*h^2*nu + 24*Rx0*D^2*R0y*h^2 + 16*D^3*R0y*h*nu^3 - 56*D^3*R0y*h*nu^2 + 48*D^3*R0y*h + 8*K0y*Rx0*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 24*Rx0*D^3*h*nu^2 + 48*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*(- 8*D^2*nu^2 + 4*Rx0*h*D*nu + 8*D^2 + 4*Rx0*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (4*D*nu)/(2*D + R0y*h) - (4*D*nu)/(2*D + Rx0*h) - (2*(8*D^2*nu + 2*D*R0y*h*nu + 2*D*Rx0*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) - (8*(- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (2*Kx0*R0y*Rx0^2*h^6 + 4*Kx0*D*Rx0^2*h^5 + 8*Kx0*R0y*D*Rx0*h^5 - 8*Kx0*D^2*Rx0*h^4*nu^2 + 16*Kx0*D^2*Rx0*h^4 - 12*R0y*D^2*Rx0*h^2*nu^2 + 24*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 + 16*D^3*Rx0*h*nu^3 - 56*D^3*Rx0*h*nu^2 + 48*D^3*Rx0*h + 8*Kx0*R0y*D^2*h^4 - 16*Kx0*D^3*h^3*nu^2 + 16*Kx0*D^3*h^3 - 24*R0y*D^3*h*nu^2 + 48*R0y*D^3*h*nu + 48*R0y*D^3*h + 16*D^4*nu^4 - 112*D^4*nu^2 + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 20) ;

            u10_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) - (8*(4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (2*(8*D^2*nu + 8*D^2 + 4*D*R0y*h + 4*D*R0y*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) + (32*D^4*nu - 96*D^4 + 96*D^4*nu^2 - 32*D^4*nu^3 - 48*D^3*R0y*h - 48*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 24*D^2*R0y*Rx0*h^2 + 48*D^3*R0y*h*nu^2 - 16*D^3*R0y*h*nu^3 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (32*D^4*nu + 64*D^4 - 96*D^4*nu^2 - 32*D^4*nu^3 + 32*D^4*nu^4 + 32*D^3*R0y*h + 32*D^3*Rx0*h + 64*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu + 16*D^2*R0y*Rx0*h^2 - 32*D^3*R0y*h*nu^2 - 16*D^3*Rx0*h*nu^2 - 16*D^2*R0y*Rx0*h^2*nu^2 + 32*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + (32*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - 8) ;

            u20_c = ((Rx0*D*R0y^2*h^3 + 2*D^2*R0y^2*h^2 + 4*Rx0*D^2*R0y*h^2 - 4*D^3*R0y*h*nu^2 + 8*D^3*R0y*h + 4*Rx0*D^3*h - 8*D^4*nu^2 + 8*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (2*(4*nu*D^2 + 2*R0y*h*nu*D))/((2*D + R0y*h)*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + (32*D^4*nu - 16*D^4*nu^2 - 32*D^4*nu^3 + 16*D^4*nu^4 + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 8*D^3*R0y*h*nu^2 - 8*D^3*Rx0*h*nu^2 - 4*D^2*R0y*Rx0*h^2*nu^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 1) ;

            u01_c = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (8*(4*D^2*nu^2 - 4*D^2 - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (2*(8*D^2*nu + 8*D^2 + 4*D*Rx0*h + 4*D*Rx0*h*nu))/((2*D + R0y*h)*(2*D + Rx0*h)) + (32*D^4*nu - 96*D^4 + 96*D^4*nu^2 - 32*D^4*nu^3 - 48*D^3*R0y*h - 48*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 24*D^2*R0y*Rx0*h^2 + 48*D^3*Rx0*h*nu^2 - 16*D^3*Rx0*h*nu^3 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (32*D^4*nu + 64*D^4 - 96*D^4*nu^2 - 32*D^4*nu^3 + 32*D^4*nu^4 + 32*D^3*R0y*h + 32*D^3*Rx0*h + 16*D^3*R0y*h*nu + 64*D^3*Rx0*h*nu + 16*D^2*R0y*Rx0*h^2 - 16*D^3*R0y*h*nu^2 - 32*D^3*Rx0*h*nu^2 - 16*D^2*R0y*Rx0*h^2*nu^2 + 32*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + (32*D*Rx0*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - 8) ;

            u02_c = ((R0y*D*Rx0^2*h^3 + 2*D^2*Rx0^2*h^2 + 4*R0y*D^2*Rx0*h^2 - 4*D^3*Rx0*h*nu^2 + 8*D^3*Rx0*h + 4*R0y*D^3*h - 8*D^4*nu^2 + 8*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (2*(4*nu*D^2 + 2*Rx0*h*nu*D))/((2*D + R0y*h)*(2*D + Rx0*h)) - (4*D*nu)/(2*D + R0y*h) + (32*D^4*nu - 16*D^4*nu^2 - 32*D^4*nu^3 + 16*D^4*nu^4 + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 8*D^3*R0y*h*nu^2 - 8*D^3*Rx0*h*nu^2 - 4*D^2*R0y*Rx0*h^2*nu^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 1) ;

            u11_c = (2 - (2*(2*D - Rx0*h))/(2*D + Rx0*h) - (2*(12*D^2 + 2*D*R0y*h + 2*D*Rx0*h - R0y*Rx0*h^2))/((2*D + R0y*h)*(2*D + Rx0*h)) - (32*D^4*nu - 64*D^4 + 64*D^4*nu^2 - 32*D^4*nu^3 - 32*D^3*R0y*h - 32*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 16*D^2*R0y*Rx0*h^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (32*D^4*nu - 64*D^4 + 64*D^4*nu^2 - 32*D^4*nu^3 - 32*D^3*R0y*h - 32*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 16*D^2*R0y*Rx0*h^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;

            out = [u00_c,u10_c,u20_c,u01_c,u02_c,u11_c];
            
        case "D01"

            u01_c  = ((4*D^2*nu^2 - 4*D^2 - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) - (4*D*nu)/(2*D + R0y*h) + (2*K0y*Rx0*R0y^2*h^6 + 4*K0y*D*R0y^2*h^5 + 8*K0y*Rx0*D*R0y*h^5 - 8*K0y*D^2*R0y*h^4*nu^2 + 16*K0y*D^2*R0y*h^4 - 14*Rx0*D^2*R0y*h^2*nu^2 + 28*Rx0*D^2*R0y*h^2*nu + 24*Rx0*D^2*R0y*h^2 - 20*D^3*R0y*h*nu^2 + 40*D^3*R0y*h*nu + 48*D^3*R0y*h + 8*K0y*Rx0*D^2*h^4 - 16*K0y*D^3*h^3*nu^2 + 16*K0y*D^3*h^3 - 28*Rx0*D^3*h*nu^2 + 56*Rx0*D^3*h*nu + 48*Rx0*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*Rx0*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20) ;

            u11_c  = ((8*(2*D - R0y*h))/(2*D + R0y*h) + (32*D^4*nu - 96*D^4 + 96*D^4*nu^2 - 32*D^4*nu^3 - 48*D^3*R0y*h - 48*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 24*D^2*R0y*Rx0*h^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

            u21_c  = ((Rx0*D*R0y^2*h^3 + 2*D^2*R0y^2*h^2 + 4*Rx0*D^2*R0y*h^2 - 4*D^3*R0y*h*nu^2 + 8*D^3*R0y*h + 4*Rx0*D^3*h - 8*D^4*nu^2 + 8*D^4)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 1) ;

            u00_c  = ((- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (2*(- 8*D^2*nu^2 + 4*Rx0*h*D*nu + 8*D^2 + 4*Rx0*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (16*D*nu)/(2*D + R0y*h) - (32*D^4*nu + 32*D^4 - 48*D^4*nu^2 - 32*D^4*nu^3 + 16*D^4*nu^4 + 16*D^3*R0y*h + 16*D^3*Rx0*h + 16*D^3*R0y*h*nu + 32*D^3*Rx0*h*nu + 8*D^2*R0y*Rx0*h^2 - 24*D^3*R0y*h*nu^2 + 8*D^3*R0y*h*nu^3 - 16*D^3*Rx0*h*nu^2 - 8*D^2*R0y*Rx0*h^2*nu^2 + 16*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

            u02_c  = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) + (16*D*nu)/(2*D + R0y*h) - (64*D^4*nu + 32*D^4 - 64*D^4*nu^2 - 64*D^4*nu^3 + 32*D^4*nu^4 + 16*D^3*R0y*h + 16*D^3*Rx0*h + 32*D^3*R0y*h*nu + 32*D^3*Rx0*h*nu + 8*D^2*R0y*Rx0*h^2 - 16*D^3*R0y*h*nu^2 - 16*D^3*Rx0*h*nu^2 - 8*D^2*R0y*Rx0*h^2*nu^2 + 16*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

            u03_c  = ((16*D^4*nu - 8*D^4*nu^2 - 16*D^4*nu^3 + 8*D^4*nu^4 + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 4*D^3*R0y*h*nu^2 - 4*D^3*Rx0*h*nu^2 - 2*D^2*R0y*Rx0*h^2*nu^2 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (4*D*nu)/(2*D + R0y*h) + 1) ;

            u12_c  = (2 - (16*D^4*nu - 32*D^4 + 32*D^4*nu^2 - 16*D^4*nu^3 - 16*D^3*R0y*h - 16*D^3*Rx0*h + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 8*D^2*R0y*Rx0*h^2 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;

            u10_c  = ((2*(4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (16*D^4*nu - 32*D^4 + 32*D^4*nu^2 - 16*D^4*nu^3 - 16*D^3*R0y*h - 16*D^3*Rx0*h + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 8*D^2*R0y*Rx0*h^2 + 16*D^3*R0y*h*nu^2 - 8*D^3*R0y*h*nu^3 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + R0y*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (4*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 2) ;

            out = [u01_c,u11_c,u21_c,u00_c,u02_c,u03_c,u12_c,u10_c];
        
        case "D02"
            u02_c = ((2*K0y*R0y*h^4 + 4*K0y*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + R0y*h)) - (8*D*nu)/(2*D + R0y*h) - (8*(4*D + 4*D*nu))/(2*D + R0y*h) + 20) ;

            u12_c = ((8*(2*D - R0y*h))/(2*D + R0y*h) + (8*D^2*nu - 24*D^2)/(D*(2*D + R0y*h)) - 8) ;

            u22_c = ((2*D^2 + R0y*h*D)/(D*(2*D + R0y*h)) + 1) ;

            u01_c = ((2*(4*D + 4*D*nu))/(2*D + R0y*h) - (- 8*D^2*nu^2 + 16*D^2*nu + 8*D^2)/(D*(2*D + R0y*h)) + (16*D*nu)/(2*D + R0y*h) - 8) ;
            u03_c = u01_c;

            u04_c = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + R0y*h)) - (4*D*nu)/(2*D + R0y*h) + 1) ;
            u00_c = u04_c;

            u13_c = (2 - (4*D^2*nu - 8*D^2)/(D*(2*D + R0y*h)) - (2*(2*D - R0y*h))/(2*D + R0y*h)) ;
            u11_c = u13_c;

            out = [u02_c,u12_c,u22_c,u01_c,u03_c,u04_c,u00_c,u13_c,u11_c];

        case "D10"
            u10_c = ((4*D^2*nu^2 - 4*D^2 + 2*D*R0y*h - 2*D*Rx0*h + R0y*Rx0*h^2)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) - (4*D*nu)/(2*D + Rx0*h) + (2*Kx0*R0y*Rx0^2*h^6 + 4*Kx0*D*Rx0^2*h^5 + 8*Kx0*R0y*D*Rx0*h^5 - 8*Kx0*D^2*Rx0*h^4*nu^2 + 16*Kx0*D^2*Rx0*h^4 - 14*R0y*D^2*Rx0*h^2*nu^2 + 28*R0y*D^2*Rx0*h^2*nu + 24*R0y*D^2*Rx0*h^2 - 20*D^3*Rx0*h*nu^2 + 40*D^3*Rx0*h*nu + 48*D^3*Rx0*h + 8*Kx0*R0y*D^2*h^4 - 16*Kx0*D^3*h^3*nu^2 + 16*Kx0*D^3*h^3 - 28*R0y*D^3*h*nu^2 + 56*R0y*D^3*h*nu + 48*R0y*D^3*h + 40*D^4*nu^4 - 80*D^4*nu^3 - 136*D^4*nu^2 + 80*D^4*nu + 96*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (8*D*R0y*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 20) ;

            u20_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) + (16*D*nu)/(2*D + Rx0*h) - (64*D^4*nu + 32*D^4 - 64*D^4*nu^2 - 64*D^4*nu^3 + 32*D^4*nu^4 + 16*D^3*R0y*h + 16*D^3*Rx0*h + 32*D^3*R0y*h*nu + 32*D^3*Rx0*h*nu + 8*D^2*R0y*Rx0*h^2 - 16*D^3*R0y*h*nu^2 - 16*D^3*Rx0*h*nu^2 - 8*D^2*R0y*Rx0*h^2*nu^2 + 16*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

            u30_c = ((16*D^4*nu - 8*D^4*nu^2 - 16*D^4*nu^3 + 8*D^4*nu^4 + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 4*D^3*R0y*h*nu^2 - 4*D^3*Rx0*h*nu^2 - 2*D^2*R0y*Rx0*h^2*nu^2 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (4*D*nu)/(2*D + Rx0*h) + 1) ;

            u00_c = ((2*(- 8*D^2*nu^2 + 4*R0y*h*D*nu + 8*D^2 + 4*R0y*h*D))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (- 8*D^2*nu^2 + 4*Rx0*h*D*nu + 8*D^2 + 4*Rx0*h*D)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + (16*D*nu)/(2*D + Rx0*h) - (32*D^4*nu + 32*D^4 - 48*D^4*nu^2 - 32*D^4*nu^3 + 16*D^4*nu^4 + 16*D^3*R0y*h + 16*D^3*Rx0*h + 32*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu + 8*D^2*R0y*Rx0*h^2 - 16*D^3*R0y*h*nu^2 - 24*D^3*Rx0*h*nu^2 + 8*D^3*Rx0*h*nu^3 - 8*D^2*R0y*Rx0*h^2*nu^2 + 16*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

            u11_c = ((8*(2*D - Rx0*h))/(2*D + Rx0*h) + (32*D^4*nu - 96*D^4 + 96*D^4*nu^2 - 32*D^4*nu^3 - 48*D^3*R0y*h - 48*D^3*Rx0*h + 16*D^3*R0y*h*nu + 16*D^3*Rx0*h*nu - 24*D^2*R0y*Rx0*h^2 + 8*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - 8) ;

            u12_c = ((R0y*D*Rx0^2*h^3 + 2*D^2*Rx0^2*h^2 + 4*R0y*D^2*Rx0*h^2 - 4*D^3*Rx0*h*nu^2 + 8*D^3*Rx0*h + 4*R0y*D^3*h - 8*D^4*nu^2 + 8*D^4)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) + 1) ;

            u21_c = (2 - (16*D^4*nu - 32*D^4 + 32*D^4*nu^2 - 16*D^4*nu^3 - 16*D^3*R0y*h - 16*D^3*Rx0*h + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 8*D^2*R0y*Rx0*h^2 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (2*(2*D - Rx0*h))/(2*D + Rx0*h)) ;

            u01_c = ((2*(4*D^2*nu^2 - 4*D^2 - 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2))/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) - (16*D^4*nu - 32*D^4 + 32*D^4*nu^2 - 16*D^4*nu^3 - 16*D^3*R0y*h - 16*D^3*Rx0*h + 8*D^3*R0y*h*nu + 8*D^3*Rx0*h*nu - 8*D^2*R0y*Rx0*h^2 + 16*D^3*Rx0*h*nu^2 - 8*D^3*Rx0*h*nu^3 + 4*D^2*R0y*Rx0*h^2*nu)/(D*(2*D + Rx0*h)*(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2)) - (4*D*Rx0*h*nu)/(4*D^2 - 4*D^2*nu^2 + 2*D*R0y*h + 2*D*Rx0*h + R0y*Rx0*h^2) + 2) ;

            out = [u10_c, u20_c, u30_c, u00_c, u11_c, u12_c, u21_c, u01_c];

        case "D11"
            out = [(20 - (2*D - Rx0*h)/(2*D + Rx0*h) - (2*D - R0y*h)/(2*D + R0y*h)), - 8, 1, ...
                ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8), ((4*D + 4*D*nu)/(2*D + R0y*h) - 8), - 8, ...
                1, 2, (2 - (2*D*nu)/(2*D + Rx0*h)), (2 - (2*D*nu)/(2*D + Rx0*h) - (2*D*nu)/(2*D + R0y*h)),...
                (2 - (2*D*nu)/(2*D + R0y*h))];

        case "D12"
            out = [(20 - (2*D - R0y*h)/(2*D + R0y*h)), - 8, 1, - 8, 1, ((4*D + 4*D*nu)/(2*D + R0y*h) - 8), - 8, 1, 2, 2, (2 - (2*D*nu)/(2*D + R0y*h)), (2 - (2*D*nu)/(2*D + R0y*h))];

        case "D20"
            u20_c = ((2*Kx0*Rx0*h^4 + 4*Kx0*D*h^3 - 12*D^2*nu^2 + 24*D^2*nu + 24*D^2)/(D*(2*D + Rx0*h)) - (8*D*nu)/(2*D + Rx0*h) - (8*(4*D + 4*D*nu))/(2*D + Rx0*h) + 20) ;

            u21_c = ((8*(2*D - Rx0*h))/(2*D + Rx0*h) + (8*D^2*nu - 24*D^2)/(D*(2*D + Rx0*h)) - 8) ;

            u22_c = ((2*D^2 + Rx0*h*D)/(D*(2*D + Rx0*h)) + 1) ;

            u10_c = ((2*(4*D + 4*D*nu))/(2*D + Rx0*h) - (- 8*D^2*nu^2 + 16*D^2*nu + 8*D^2)/(D*(2*D + Rx0*h)) + (16*D*nu)/(2*D + Rx0*h) - 8) ;
            u30_c = u10_c;

            u40_c = ((- 2*D^2*nu^2 + 4*D^2*nu)/(D*(2*D + Rx0*h)) - (4*D*nu)/(2*D + Rx0*h) + 1) ;
            u00_c = u40_c;

            u31_c = (2 - (4*D^2*nu - 8*D^2)/(D*(2*D + Rx0*h)) - (2*(2*D - Rx0*h))/(2*D + Rx0*h)) ;
            u11_c = u31_c;

            out = [u20_c,u21_c,u22_c,u10_c,u30_c,u40_c,u00_c,u31_c,u11_c];

        case "D21"
            out = [(20 - (2*D - Rx0*h)/(2*D + Rx0*h)), - 8, 1, + ((4*D + 4*D*nu)/(2*D + Rx0*h) - 8),...
           - 8, - 8, 1, 1, 2, (2 - (2*D*nu)/(2*D + Rx0*h)), (2 - (2*D*nu)/(2*D + Rx0*h)), 2];

        case "D22"
            out = [1, 2, -8, 2, 1, -8, 20, -8, 1, 2, -8, 2, 1];
    end
    varargout = num2cell(out);
end

