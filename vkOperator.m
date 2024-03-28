function L = vkOperator(F1,F2,Dx,Dy,Dxx,Dyy)

% Dxx = Dx*Dx ;
% Dyy = Dy*Dy ;
Dxy = Dx*Dy ;

L = (Dxx*F1).*(Dyy*F2) + (Dyy*F1).*(Dxx*F2) - 2*(Dxy*F1).*(Dxy*F2) ;

end
