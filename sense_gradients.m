% takes in a cell (x_cell,y_cell) position and lattice of chemoattractant
% concentrations (ca, on lattice x, y) and outputs an angle, theta, of
% gradient signal sensed

% based on Louise Dyson D.Phil project 
% extended by L.J.Schumacher


function [filopodia,move,theta,num_better,deltaC] = sense_gradients(sampledDirections,...
    x_cell,y_cell,ca,x,y,eatWidth,filolength,numFilopodia,fil,sensingAccuracy)


%% Integrate the chemoattractant in the present area %%%
multiplier = 1/(eatWidth^2*2*pi)*exp(-(x-x_cell).^2/2/eatWidth^2)*exp(-(y'-y_cell).^2/2/eatWidth^2);

%% Trapezoid sum computed with vector-matrix multiply. %
intgrand = ca.*multiplier;
intgrand = diff(x')*(intgrand(1:end-1,:) + intgrand(2:end,:))./2; 
present_area = (intgrand(1:end-1) + intgrand(2:end))*diff(y)./2; % include the lattice spacing to get correct integral result -- LJS

%% Integrate the chemoattractant in the area of the filopodia %%%
new_area = 0;
num_better = 0;
chosen_theta_idx = NaN;
filopodia = NaN(numFilopodia,2); % store the x-y coordinates of all the cell's filopodia -- LJS
for filo_ctr=1:numFilopodia %loops through the filopodia and keeps track of best direction -- LJS
    if isempty(fil)
        %% the cell extends a filopodia in the sampledDirections direction %%%
        x_fil = x_cell + cos(sampledDirections(filo_ctr))*filolength;  % x coordinate of the filopodia
        y_fil = y_cell + sin(sampledDirections(filo_ctr))*filolength;  % y coordinate of the filopodia
        filopodia(filo_ctr,:) = [x_fil,y_fil];
    else
        x_fil = fil(filo_ctr,1);
        y_fil = fil(filo_ctr,2);
        filopodia(filo_ctr,:) = fil(filo_ctr,:);
        sampledDirections(filo_ctr) = atan2(y_fil - y_cell, x_fil - x_cell);
    end
    old_area = new_area;
    multiplier = 1/(eatWidth^2*2*pi)*exp(-(x-x_fil).^2/2/eatWidth^2)*exp(-(y'-y_fil).^2/2/eatWidth^2);

%%    Trapezoid sum computed with vector-matrix multiply. %
    intgrand = ca.*multiplier;
    intgrand = diff(x')*(intgrand(1:end-1,:) + intgrand(2:end,:))./2; 
    new_area = (intgrand(1:end-1) + intgrand(2:end))*diff(y)./2; % include the lattice spacing to get correct intergral result -- LJS

    if new_area>present_area
        num_better = num_better+1;
    end
    if new_area>=old_area
        chosen_theta_idx = filo_ctr;
    else
        new_area = old_area;
    end
end

%% If the present area is better then stay put, else choose to move in the theta direction%%%
deltaC = (new_area - present_area);
caDiff = deltaC/present_area;    % if caDiff > sensingAccuracy/sqrt(ca) then the cell moves in sensed direction -- LJS
if caDiff < sensingAccuracy/sqrt(present_area) % then present_area > new_area and the cell doesn't try to move
    move = 0;
    theta = (rand()*2 - 1)*pi; % direction of movement is random -- LJS
else   % else the cell does try to move
    move=1;
    if ~isempty(sampledDirections)
        theta = sampledDirections(chosen_theta_idx);
    else
        theta = (rand()*2 - 1)*pi;% direction of movement is random -- LJS
    end
end
