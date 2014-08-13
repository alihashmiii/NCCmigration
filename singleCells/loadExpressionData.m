% load gene expression data of single cells into array of genes by cell
% number for each cell group (16hr, 24h trailblazers, 24h quartiles 1-4)

raw = importdata('16 and 24 hr trailblazers with 24hr quartiles individual cell log2Ex for PARTEK 063014.csv');
genes = raw.textdata(2:end,1);
trailblazers16h = raw.data(:,2:73);
trailblazers24h = raw.data(:,74:149);
quartile1 = raw.data(:,150:192);
quartile2 = raw.data(:,193:233);
quartile3 = raw.data(:,234:277);
quartile4 = raw.data(:,278:318);

clear raw