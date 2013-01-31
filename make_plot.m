%% A function to plot the data in cells, cells_follow and ca

function [] = make_plot(cells,cells_follow,xlat,ylat,ca,filopodia,attach,cell_radius,edge,barrier,experiment)

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
if cell_radius < 1
    filopodia = cells';
    cell_radius = 5;
end

if (isempty(filopodia)==1)
    filopodia = cells';
end

t = 0:0.1:2*pi;
for cell_ind = 1:length(cells(1,:))
    if cells_follow(cell_ind)==1
        if attach(cell_ind)==0
            fill(cell_radius*cos(t)+cells(1,cell_ind),cell_radius*sin(t)+cells(2,cell_ind),'r');
            plot([cells(1,cell_ind),filopodia(cell_ind,1)],[cells(2,cell_ind),filopodia(cell_ind,2)],'r','LineWidth',2)
        else
            fill(cell_radius*cos(t)+cells(1,cell_ind),cell_radius*sin(t)+cells(2,cell_ind),'w');
            plot([cells(1,cell_ind),filopodia(cell_ind,1)],[cells(2,cell_ind),filopodia(cell_ind,2)],'w','LineWidth',2)
        end
    else
        fill(cell_radius*cos(t)+cells(1,cell_ind),cell_radius*sin(t)+cells(2,cell_ind),'y');
        hold on
        plot([cells(1,cell_ind),filopodia(cell_ind,1)],[cells(2,cell_ind),filopodia(cell_ind,2)],'y','LineWidth',2)
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