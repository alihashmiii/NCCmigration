function neighbours = neighbourRelationships(cells)
% neighbour relationships 
% input cells should be of shape [x1, y1; x2, y2; ...; xN, yN] where N is
% number of cells
if size(cells,1) < size(cells,2)
    cells = double(cells');
else
    cells = double(cells);
end
numCells = size(cells,1);
neighbourCutoff = 84.34; % in mu, if cells are further apart don't count them as neighbours
areaCutoff = pi*(neighbourCutoff/2).^2; % important to cut out voronoi areas e.g. at the edges

DT = delaunayTriangulation(cells);
cellEdges = edges(DT); %edge indices
distXY = cells(cellEdges(:,1),:) - cells(cellEdges(:,2),:);
edgeLengths = sqrt(distXY(:,1).^2 + distXY(:,2).^2);
neighbourIndcs = edgeLengths<=neighbourCutoff;
neighbourEdges = cellEdges(neighbourIndcs,:);

cellRadius = 7.5;
neighbours.distancesBinEdges = 0:cellRadius:neighbourCutoff;
neighbours.distances = histc(edgeLengths(neighbourIndcs),neighbours.distancesBinEdges);
neighbours.numbers = hist(hist(neighbourEdges(:),numCells),1:10);

[voronoiVertices, voronoiRegions] = voronoiDiagram(DT);

% loop through remaining cells to calculate their neighbourhood area
voronoiAreas = NaN(numCells,1);
excludedCells = false(numCells,1);
for cellCtr = 1:numCells
    cellRegionVertices = voronoiVertices(voronoiRegions{cellCtr},:);
    voronoiAreas(cellCtr) = polyarea(cellRegionVertices(:,1), ...
        cellRegionVertices(:,2)); % area for polygons with nodes at inf calculates to NaN
    
    % calculate node-vertex distances (for excluding edge-distorted
    % regions)
    distNV = sqrt((cellRegionVertices(:,1) - cells(cellCtr,1)).^2 ...
        + (cellRegionVertices(:,2) - cells(cellCtr,2)).^2);
    if any(distNV>neighbourCutoff)
        excludedCells(cellCtr) = true;
    end
end
% exclude voronoi regions that contain delaunay-edges longer than the
% cut-off (~neighbourIndcs), or that have cells on the convex hull (these have vertices at infity
% and hence NaN area), or that have node-vertex distances greater than the
% neighbourCutoff (excludedCells)
excludedCellIndcs = unique([unique(cellEdges(~neighbourIndcs,:)); ...
    convexHull(DT); find(excludedCells)]);
% also exclude any voronoi polygons with too big an area
neighbours.areasBinEdges = linspace(0,areaCutoff,21);
neighbours.areas = histc(...
    voronoiAreas(voronoiAreas<=areaCutoff&...
    ~ismember(1:numCells,excludedCellIndcs)')...
    ,neighbours.areasBinEdges);
