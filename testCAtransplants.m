% test VEGF transplant simulations
in.time = 18;
% run solver with and without cells
for insertCells = 0:1
    in.insertCells = insertCells;
    for experiment = 11:14
        if ~insertCells
            in.numCellsInitial = 1;
        else
            in.numCellsInitial = 6;
        end
        in.experiment = experiment;
        CA6(in,experiment)
    end
end