function D4 = beam_biharmonic_build(M,L,bc)

%---- this function creates the 1D biharmonic operator over a domain of length
%---- L, using M grid intervals, and using 3 different kinds of boundary
%---- conditions: 1: Free Free, 2: SimS SimS, 3: Clamp Clamp


h = L/M ;
v = ones(M+1,1) ;


if bc == 1
    
    D2 = spdiags([v,-2*v,v],-1:1,M+1,M+1) ;
    D4 = D2*D2 ;

    D4(1,1) = 1 ; D4(1,2) = -2;
    D4(end,end) = 1 ; D4(end,end-1) = -2;
    D4(2,1) = -2 ; D4(2,2) = 5 ; D4(2,3) = -4 ;
    D4(end-1,end) = -2 ; D4(end-1,end-1) = 5 ; D4(end-1,end-2) = -4 ;


elseif bc == 2

    D2 = spdiags([v,-2*v,v],-1:1,M-1,M-1) ;
    D4 = D2*D2 ;


elseif bc == 3

    D4 = spdiags([v,-4*v,6*v,-4*v,v],-2:2,M-3,M-3) ;
   % D4(1) = 6 ; D4(end) = 6 ;

end

D4        = D4/h^4 ;

end


