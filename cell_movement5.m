% Louise Dyson D.Phil project second CA program, 22/10/09
% takes in a cell (x_cell,y_cell) position and lattice of chemoattractant
% concentrations (ca, on lattice x, y) and outputs an angle, theta, of
% movement for the cell and a parameter, d, that is related to the radius
% of sensing of ca.

function [filopodia,caDiff,move,theta,num_better] = cell_movement5(theta,x_cell,y_cell,ca,x,y,eatWidth,filolength,numFilopodia,fil)


%% Integrate the chemoattractant in the present area %%%
multiplier = 1/(eatWidth^2*2*pi)*exp(-(x-x_cell).^2/2/eatWidth^2)*exp(-(y'-y_cell).^2/2/eatWidth^2);

%% Trapezoid sum computed with vector-matrix multiply. %
intgrand = ca.*multiplier;
intgrand = sum((intgrand(1:end-1,:) + intgrand(2:end,:))./2,1);
present_area = sum((intgrand(1:end-1) + intgrand(2:end))./2);

%% Integrate the chemoattractant in the area of the filopodia %%%
new_area = 0;
num_better = 0;
chosen_theta_idx = NaN;
filopodia = NaN(numFilopodia,2); % store the x-y coordinates of all the cell's filopodia -- LJS
for filo_ctr=1:numFilopodia %loops through the filopodia and keeps track of best direction -- LJS
    if isempty(fil)
        %% the cell extends a filopodia in the theta direction %%%
        x_fil = x_cell + cos(theta(filo_ctr))*filolength;  % x coordinate of the filopodia
        y_fil = y_cell + sin(theta(filo_ctr))*filolength;  % y coordinate of the filopodia
        filopodia(filo_ctr,:) = [x_fil,y_fil];
    else
        x_fil = fil(filo_ctr,1);
        y_fil = fil(filo_ctr,2);
        filopodia(filo_ctr,:) = fil(filo_ctr,:);
    end
    old_area = new_area;
    multiplier = 1/(eatWidth^2*2*pi)*exp(-(x-x_fil).^2/2/eatWidth^2)*exp(-(y'-y_fil).^2/2/eatWidth^2);

%%    Trapezoid sum computed with vector-matrix multiply. %
    intgrand = ca.*multiplier;
    intgrand = sum((intgrand(1:end-1,:) + intgrand(2:end,:))./2,1);
    new_area = sum((intgrand(1:end-1) + intgrand(2:end))./2);
    if new_area>present_area
        num_better = num_better+1;
    end
    if new_area>=old_area
        chosen_theta_idx = filo_ctr;
    else
        new_area = old_area;
    end
end

%% If the present area is better then stay put, else move in the theta direction%%%
caDiff = (present_area - new_area)/present_area;    % energy difference. So if caDiff<0 then the cell definitely moves
if caDiff>0 % then present_area>new_area and the cell doesn't try to move
    move=0;
    theta = NaN; % direction of movement is undefined -- LJS
else   % else the cell does try to move
    move=1; 
    theta = theta(chosen_theta_idx);
end
