%% A function to plot the data in cells, cellsFollow and ca

function [] = make_plot(cells,cellsFollow,xlat,ylat,ca,filopodia,num_filopodia,attach,cellRadius,edge,barrier,experiment)

if experiment~=6
    if isempty(ca)==1
        if isempty(filopodia)==1    % for make_comp_movie.m
            filopodia = cells';
            whitebg('white')
            set(gca,'color',[0.4,0.4,0.4]);
        else
            whitebg('blue')
        end
    elseif edge==0
        whitebg('white')
        contourf(xlat,ylat,ca',10,'EdgeColor','none')
        hold on
    else
        whitebg('white')
        contourf(xlat,ylat,ca',10,'EdgeColor','none')
        hold on
    end
end
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
        filoNum = num_filopodia(2);
        if attach(cellidx)==0
            cellColor = 'r';
        else
            cellColor = 'w';
        end
    else
        filoNum = num_filopodia(1);
        cellColor = 'y';
    end
    
    % draw the cell
    fill(cellRadius*cos(t)+cells(1,cellidx),cellRadius*sin(t)+cells(2,cellidx),cellColor);
    
    % draw the filopodia
    for filoidx = 1:filoNum
        plot([cells(1,cellidx),filopodia(cellidx,filoidx,1)],[cells(2,cellidx),filopodia(cellidx,filoidx,2)],cellColor,'LineWidth',1)
    end
end

if (experiment==4)||(experiment==5)
    hold on
    y = 0:120;
    plot(barrier.*ones(length(y)),y)
elseif (experiment==6)
    hold on
    y = 0:120;
    plot(max(xlat).*ones(length(y)),y)
end
hold off
axis equal
colorbar, caxis([0 max(max(ca))])

ylim([0,120])
set(gca,'YTick',[0,120])
xlim([min(xlat),1000])