% function to plot cells in domain with chemoattractant gradient

function [] = make_plot(cells,cellsFollow,xlat,ylat,time,ca,filopodia,numFilopodia,attach,cellRadius,filolength,sensingAccuracy,showColorbar,caCmap,quickMode,param,dan)
experiment = param.experiment;
transplantTime = param.transplantTime;
transplantXLocation = param.transplantXLocation;
domainHeight = param.domainHeight;

leadColor = [1 1 0.5];%[247 252 185]/255;
followColor = [111 192 113]/255;
lostColor = [0.5 0.5 0.5];

if ~quickMode, whitebg('white'), end
%  plot the chemoattractant
contourf(xlat,ylat,ca',100,'EdgeColor','none')
hold on

% % calculate and plot the CA gradient in regions where it can be sensed
% [dcadx, dcady] = gradient(ca',xlat,ylat);
% indicesSensible = sqrt(dcadx.^2 + dcady.^2)*filolength./sqrt(ca')>=sensingAccuracy;
% if any(~indicesSensible(:))
%     % plot a shaded region of (non-)sensible region of CA
%     sensAccRegion = pcolor(xlat,ylat,zeros(length(ylat),length(xlat)));
%     set(sensAccRegion,'EdgeColor','none','FaceColor',[0.75 0.75 0.75],...
%         'AlphaData',0.5*double(~indicesSensible),'FaceAlpha','interp','AlphaDataMapping','none');
% %     % plot markers on the lattice points in the non-sensible region
% %     [yIndcsSnsbl, xIndcsSnsbl] = find(~indicesSensible);
% %     scatter(xlat(xIndcsSnsbl),ylat(yIndcsSnsbl),20,[0.75 0.75 0.75],'x');
%     % draw a contour around the non-sensible region
%     sensAccContour = contour(xlat,ylat,indicesSensible,1,'EdgeColor',[0.75 0.75 0.75], 'LineWidth', 2);
% end

% if doing a VEGF transplant simulation, plot regions of increased CA production
if (experiment==11||experiment==12||experiment==13||experiment==14)&&time>=transplantTime
    if (experiment==12||experiment==13)
        xindcs =  (xlat>=transplantXLocation)&(xlat<=(transplantXLocation + 1/8*max(xlat)));
        yindcs = (ylat<=domainHeight/2);
    elseif experiment==11 % increased chemoattractant at edge of domain, close to entrance
        xindcs =  (xlat>=transplantXLocation)&(xlat<=(transplantXLocation + 1/8*max(xlat)));
        yindcs = (ylat<=domainHeight/20);
    elseif experiment==14 % increased chemoattractant at far end of domain, full width
        xindcs =  (xlat>=transplantXLocation)&(xlat<=max(xlat));
        yindcs = true(size(ylat));
    end
    contour(xlat,ylat,single(yindcs)*single(xindcs)',1,'EdgeColor',[1 0.2 0.2],'LineWidth',2)
end

if (experiment==40)||(experiment==41)||(experiment==42)||(experiment==43)||experiment==44 % for DAN simulations, show slow zone
    xmax = max(xlat)/3;
    ymin = min(ylat);
    ymax = max(ylat);
    if experiment==42||experiment==44 % determine slow-down for dynamic dan for color and hatch spacing
        tPeakSlowdown = 12;
        minSlowdown = 0.5;
        slowDown = max(minSlowdown,(tPeakSlowdown - abs(time - tPeakSlowdown))/tPeakSlowdown);
    elseif experiment==41 % dilutiuon of DAN through growth
        slowDown = param.initialDomainLength/max(xlat);
    else
        slowDown = 1;
    end
    patchHandle = patch([0,max(xlat)/3,max(xlat)/3,0,0],[ymin,ymin,ymax,ymax,ymin],...
        'r','FaceColor','none','EdgeColor',[1 0.2 0.2],'LineWidth',2);
    hatchHandle = hatchfill2(patchHandle,'HatchSpacing',10/slowDown,'HatchLineWidth',2);
    if experiment==43||experiment==44 % show contour of unconsumed dan
        [~, contourHandle] = contour(xlat(xlat<=xmax+eps(xmax)),ylat,dan,1,...
            'EdgeColor',[0.5 0.5 0.5],'FaceColor','none','LineWidth',2);
        hatchfill2(contourHandle,'HatchSpacing',10/slowDown,'HatchLineWidth',2);
    end
    if experiment==41||experiment==42||experiment==44
        hatchHandle.Color = hatchHandle.Color*slowDown;
        patchHandle.EdgeAlpha = slowDown;
    end
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
phi = 0:0.1:2*pi; % for plotting the cell circles -- LJS

for cellidx = 1:length(cells(1,:))
    if cellsFollow(cellidx)==1
        filoNum = numFilopodia(end);
        if attach(cellidx)==0
            cellColor = lostColor; % dettached followers -- LJS
        else
            cellColor = followColor; % followers -- LJS
        end
    elseif cellsFollow(cellidx)==0
        filoNum = numFilopodia(1);
        cellColor = leadColor; % leaders
    else % continuous states
        filoNum = numFilopodia;
        if isnan(cellsFollow(cellidx))
            cellColor = lostColor;
        else
            cellColor = leadColor*(1 - cellsFollow(cellidx)) ...
                + cellsFollow(cellidx)*followColor;
        end
    end
    
    % draw the cell
    fill(cellRadius*cos(phi)+cells(1,cellidx),cellRadius*sin(phi)+cells(2,cellidx),cellColor,'EdgeColor','none');
    
    if ~quickMode
        % draw the filopodia
        for filoidx = 1:filoNum
            plot([cells(1,cellidx),filopodia(cellidx,filoidx,1)],[cells(2,cellidx),filopodia(cellidx,filoidx,2)],'color',cellColor,'LineWidth',1)
        end
    end
end

hold off
axis equal
if showColorbar
    hcb = colorbar;
    set(hcb,'ytick',[0 1])
    set(get(hcb,'title'),'string','CA')
end
set(gca,'clim',[0 1]) % if set to [0 max(max(ca))], can cause jittering of colorbar labels in movies


ylim([0,max(ylat)])
set(gca,'YTick',[0,max(ylat)])
% xlim([min(xlat),1000])
set(gca,'XTick',0:200:1000,'FontSize',12)