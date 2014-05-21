%% A function to plot the data in cells, cellsFollow and ca

function [] = make_plot(cells,cellsFollow,xlat,ylat,ca,filopodia,numFilopodia,attach,cellRadius,filolength,sensingAccuracy,showColorbar,caCmap,quickMode)

if ~quickMode, whitebg('white'), end
%  plot the chemoattractant
contourf(xlat,ylat,ca',20,'EdgeColor','none')
hold on

% calculate and plot the CA gradient in regions where it can be sensed
[dcadx, dcady] = gradient(ca',xlat,ylat);
indicesSensible = sqrt(dcadx.^2 + dcady.^2)*filolength./sqrt(ca')>=sensingAccuracy;
if any(~indicesSensible(:))
    % plot a shaded region of sensible region of CA
    sensAccRegion = pcolor(xlat,ylat,zeros(length(ylat),length(xlat)));
    set(sensAccRegion,'EdgeColor','none','FaceColor',[0.75 0.75 0.75],...
        'AlphaData',0.5*double(~indicesSensible),'FaceAlpha','interp','AlphaDataMapping','none');
    % plot markers on the lattice points in the non-sensible region
    [yIndcsSnsbl, xIndcsSnsbl] = find(~indicesSensible);
    scatter(xlat(xIndcsSnsbl),ylat(yIndcsSnsbl),20,[0.75 0.75 0.75],'x');
    % draw a contour around the non-sensible region
    sensAccContour = contour(xlat,ylat,indicesSensible,1,'EdgeColor',[0.75 0.75 0.75], 'LineWidth', 2);
else
    pause()
end

colormap(caCmap)

if cellRadius < 1
    filopodia = cells';
    cellRadius = 5;
end

if (isempty(filopodia)==1)
    filopodia = cells';
end

% draw the cells with their filopodia -- LJS
t = 0:0.1:2*pi; % for plotting the cell circles -- LJS

for cellidx = 1:length(cells(1,:))
    if cellsFollow(cellidx)==1
        filoNum = numFilopodia(2);
        if attach(cellidx)==0
            cellColor = [0.5 0.5 0.5];
        else
            cellColor = 'w';
        end
    else
        filoNum = numFilopodia(1);
        cellColor = 'y';
    end
    
    % draw the cell
    fill(cellRadius*cos(t)+cells(1,cellidx),cellRadius*sin(t)+cells(2,cellidx),cellColor);
    
    if ~quickMode
        % draw the filopodia
        for filoidx = 1:filoNum
            plot([cells(1,cellidx),filopodia(cellidx,filoidx,1)],[cells(2,cellidx),filopodia(cellidx,filoidx,2)],'color',cellColor,'LineWidth',1)
        end
    end
end

hold off
axis equal
if showColorbar, colorbar, end
set(gca,'clim',[0 1]) % if set to [0 max(max(ca))], can cause jittering of colorbar labels in movies

ylim([0,120])
set(gca,'YTick',[0,120])
xlim([min(xlat),1000])