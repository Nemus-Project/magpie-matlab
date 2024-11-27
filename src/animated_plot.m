function animated_plot(plotvec,plotScale,xax,yax,Nvec,Lvec,absPlot,plotTitle,Nribs,beamCoord,beamParams,Nlump,lumpCoord,Nouts,outMat,x_f)


Nx = Nvec(1) ;
Ny = Nvec(2) ;
Lx = Lvec(1) ;
Ly = Lvec(2) ;

plotDist = reshape(plotvec,Ny+1,Nx+1)/plotScale ;
plotC = plotDist ;
if absPlot
    plotC = abs(plotDist) ;
end

surf(xax,yax,plotDist,plotC,'LineStyle','none') ; view(2) ; axis equal ; axis tight ; colorbar; title(plotTitle,'interpreter','latex')
hold on ;
for nR = 1 : Nribs
    Nbeam = beamParams(nR,5) ;
    x_coord_beam = linspace(beamCoord(nR,1),beamCoord(nR,2),Nbeam+1) ;
    y_coord_beam = linspace(beamCoord(nR+Nribs,1),beamCoord(nR+Nribs,2),Nbeam+1) ;
    plot3(x_coord_beam*Lx,y_coord_beam*Ly,10*ones(1,Nbeam+1),'r') ;
end
for nL = 1 : Nlump
    x_coord_lump = lumpCoord(nL) ; y_coord_lump = lumpCoord(nL+Nlump) ;
    plot3( x_coord_lump*Lx, y_coord_lump*Ly,1e3,'r','marker','o','markersize',5,'markerfacecolor','w') ;
end
for nO = 1 : Nouts
    x_coord_out = outMat(nO,1) ; y_coord_out = outMat(nO,2) ;
    plot3(x_coord_out*Lx, y_coord_out*Ly,1e3,'k','marker','x','markersize',5,'markerfacecolor','b') ;
end
x_coord_in = x_f(1) ; y_coord_in = x_f(2) ;
plot3(x_coord_in*Lx, y_coord_in*Ly,1e3,'w','marker','square','markersize',5,'markerfacecolor','r') ;
hold off ;