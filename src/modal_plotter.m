function  modal_plotter(cmap,Nplots,NN,X,Y,Qplate,Lvec,Nvec,Nribs,Nlump,beamParams,beamCoord,lumpCoord)


if Nplots ~= 9 && Nplots ~= 6 && Nplots ~= 3
    warning('Number of plots set to 3. Alternatively, specify Nplots = 6, Nplots = 9') ;
end


Nx = Nvec(1) ; Ny = Nvec(2) ;
Lx = Lvec(1) ; Ly = Lvec(2) ; 

x_lump_coord = lumpCoord(1:Nlump) ;
y_lump_coord = lumpCoord(Nlump+1:end) ;
colormap(cmap) ;


%-------------------------------------------------------------------------
% PLOTS
if Nplots == 9

    for m = 1 : Nplots
        mdShape = reshape(Qplate(:,m+NN-1),[(Ny+1),(Nx+1)]) ;
        subplot(3,3,m)
        mesh(X,Y,real(mdShape),real((mdShape)),'FaceColor','texturemap') ; view(2); axis equal; axis tight; hold on ;
        for nR = 1 : Nribs
            Nbeam = beamParams(nR,5) ;
            x_coord_beam = linspace(beamCoord(nR,1),beamCoord(nR,2),Nbeam+1) ;
            y_coord_beam = linspace(beamCoord(nR+Nribs,1),beamCoord(nR+Nribs,2),Nbeam+1) ;
            plot3(x_coord_beam*Lx,y_coord_beam*Ly,10*ones(1,Nbeam+1),'r') ;
        end
        for nL = 1 : Nlump
            plot3(x_lump_coord(nL)*Lx,y_lump_coord(nL)*Ly,10,'r','marker','o','markersize',5,'markerfacecolor','w') ;
        end

        xlabel('x (m)','interpreter','latex')
        ylabel('y (m)','interpreter','latex')
        tit = sprintf('Mode %d',m+NN-1) ;
        title(tit) ;
    end
    %-------------------------------------------------------------------------


elseif Nplots == 6


    for m = 1 : Nplots
        mdShape = reshape(Qplate(:,m+NN-1),[(Ny+1),(Nx+1)]) ;
        subplot(2,3,m)
        mesh(X,Y,real(mdShape),real((mdShape)),'FaceColor','texturemap') ; view(2); axis equal; axis tight; hold on ;
        for nR = 1 : Nribs
            Nbeam = beamParams(nR,5) ;
            x_coord_beam = linspace(beamCoord(nR,1),beamCoord(nR,2),Nbeam+1) ;
            y_coord_beam = linspace(beamCoord(nR+Nribs,1),beamCoord(nR+Nribs,2),Nbeam+1) ;
            plot3(x_coord_beam*Lx,y_coord_beam*Ly,10*ones(1,Nbeam+1),'r') ;
        end
        for nL = 1 : Nlump
            plot3(x_lump_coord(nL)*Lx,y_lump_coord(nL)*Ly,10,'r','marker','o','markersize',5,'markerfacecolor','w') ;
        end

        xlabel('x (m)','interpreter','latex')
        ylabel('y (m)','interpreter','latex')
        tit = sprintf('Mode %d',m+NN-1) ;
        title(tit) ;
    end
    %-------------------------------------------------------------------------


elseif Nplots == 3


    for m = 1 : Nplots
        mdShape = reshape(Qplate(:,m+NN-1),[(Ny+1),(Nx+1)]) ;
        subplot(1,3,m)
        mesh(X,Y,real(mdShape),real((mdShape)),'FaceColor','texturemap') ; view(2); axis equal; axis tight; hold on ;
        for nR = 1 : Nribs
            Nbeam = beamParams(nR,5) ;
            x_coord_beam = linspace(beamCoord(nR,1),beamCoord(nR,2),Nbeam+1) ;
            y_coord_beam = linspace(beamCoord(nR+Nribs,1),beamCoord(nR+Nribs,2),Nbeam+1) ;
            plot3(x_coord_beam*Lx,y_coord_beam*Ly,10*ones(1,Nbeam+1),'r') ;
        end
        for nL = 1 : Nlump
            plot3(x_lump_coord(nL)*Lx,y_lump_coord(nL)*Ly,10,'r','marker','o','markersize',5,'markerfacecolor','w') ;
        end

        xlabel('x (m)','interpreter','latex')
        ylabel('y (m)','interpreter','latex')
        tit = sprintf('Mode %d',m+NN-1) ;
        title(tit) ;
    end
    %-------------------------------------------------------------------------

end

