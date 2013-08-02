% Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/03/10
% Experimental modelling for use with CA5.m

% ablates (deletes) the cells behind x_start, also dealing with the cell_followers

function out = ablate(cells,attach,x_start,direction,cellsFollow)

if strcmp(direction,'back')==1
    ablate = find(cells(1,:)<x_start);
elseif strcmp(direction,'forwards')==1
    ablate = find(cells(1,:)>x_start);
else
    error([mat2str(direction),'is not a valid input to function ablate.m'])
end
if any(cellsFollow)
    % dettach cells from ablated cells
    for i=1:length(ablate)
        attach(attach==ablate(i))=0;
        attach(attach>ablate(i))=attach(attach>ablate(i))-1;
    end
    attach(cells(1,:)<x_start)=[];      % delete attach details for ablated cells
end
cellsFollow(ablate)=[];    % delete follower details for ablated cells
cells(:,ablate)=[];         % ablate the non-front cells
out.cells = cells;
out.attach = attach;
out.cellsFollow = cellsFollow;
disp('cells ablated')