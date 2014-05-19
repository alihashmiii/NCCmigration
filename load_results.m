% script to load results and transfer from output structure fields
% to variables, so that other scripts can be run on it as if simulation had
% just been done 
% 19.05.14 -- LJS

varNames = fields(out);

for varCtr = 1:length(varNames)
    eval([ cell2mat(varNames(varCtr)) '= out.' cell2mat(varNames(varCtr)) ';']);
end