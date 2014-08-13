%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find_MI_relations.m
%
% Find_MI_1.0 - Extract mutual information relations from finite data, version 1.0.
% Copyright (C) Dec. 2004 Noam Slonim, Gurinder S. Atwal, Gasper Tkacik, and William Bialek
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%
% Introduction %
%%%%%%%%%%%%%%%%
%
% This code implements the mutual information estimation procedure used by
% Slonim et. al in the paper entitled "Information based clustering", 2005.
% Given a set of patterns, represented as the rows of the input matrix M,
% we estimate the mutual information between every pair of patterns. The
% idea behind the algorithm is based on the "direct" estimation method,
% originaly developed by Strong et. al, PRL, 1998, for the analysis of
% neural coding. This method is extended here so it can be equally applied
% to general data. See "Estimating mutual information and multi information
% in large networks", Slonim et. al, for the specific details of the
% implementation. 
% 
%%%%%%%%%
% INPUT %
%%%%%%%%%
%
% M: A matrix in which every ROW corresponds to some data pattern, e.g., the expression
% profile of a single gene under many conditions. 
%
% b_star: A scalar which sets the maximal quantization level. The default is floor(.5*sqrt(Ns)), where Ns is
% the sample size, i.e., the number of columns in M. 
%
% RandMode: A scalar flag. When its ON (i.e., not zero) the data is scrambled to remove any
% potential real correlations. In this case the resulting information
% relations should be around zero. This mode is useful to determine whether
% b_star is not too high and the results are not over estimated. The
% default is 0. 
%
% SampSizePerc: A row vector which determines the sub-sample sizes used for the extrapolation
% curve. The default is [.7 .7875 .9 1]. See the corresponding paper for
% the reasoning behind these numbers. 
% 
% Trials: A row vector which sets the number of times we estimate the information in each
% sub-sample size for the extrapolation curve. The default is [21 16 12 1]. 
% See the corresponding paper for the reasoning behind these numbers.
% Notice, that the dimension of this argument must be identical to the
% dimensions of SampSizePerc.
%
% MISS_VAL: a scalar, typically inf, which marks missing values in the
% matrix M. For every pair of rows in M, only values which are not MISS_VAL
% in BOTH rows will be used for the estimation. The default is inf. 
% 
% MinSample: A scalar which sets the minimal sample size for estimation.
% For pairs with a joint sample size which is smaller than MinSample, no
% estimation is performed. The default is 100.
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Command line examples %
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (1) [Res] = Find_MI_relations (M,[],[],[],[],[],[]);
% (2) [Res] = Find_MI_relations (M,5,[],[],[],[],[]);
% (3) [Res] = Find_MI_relations (M,5,0,[],[],[],[]);
% (4) [Res] = Find_MI_relations (M,5,0,[.7 .7875 .9 1],[],[],[]);
% (5) [Res] = Find_MI_relations (M,5,1,[.7 .7875 .9 1],[21 16 12 1],[],[]);
% (6) [Res] = Find_MI_relations (M,5,0,[.7 .7875 .9 1],[21 16 12 1],inf,[]);
% (7) [Res] = Find_MI_relations (M,5,0,[.7 .7875 .9 1],[21 16 12 1],inf,200);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%
% OUTPUT %
%%%%%%%%%%
%
% Res.I: Obtained pairwise mutual information relations for all patterns in
% M. This is an N x N matrix where N is the number of rows in M.
%
% Res.SampSize: Indicates the joint sample size for every pair (again, 
% an N x N matrix). 
%
% Res.prm: Values of the parameters used throughout the estimation. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Res] = Find_MI_relations (M,b_star,RandMode,SampSizePerc,Trials,MISS_VAL,MinSample)

Res = InitRes (M,b_star,RandMode,SampSizePerc,Trials,MISS_VAL,MinSample);

Res = QuantizeM (M,Res); 

for i=1:Res.prm.N-1    
    for j=i+1:Res.prm.N
        Res = FindPairInfo (i,j,Res);        
    end
    if mod(i,10)==0,
        fprintf ('First %d rows are done...\n',i);
    end    
end

Res = rmfield (Res,'GoodInds');
Res = rmfield (Res,'QM');    
Res.prm.RunEnd = datestr(now);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = InitRes (M,b_star,RandMode,SampSizePerc,Trials,MISS_VAL,MinSample)

Res.prm.RandMode = RandMode;
if RandMode
    Res.WARNING = 'RAND MODE ON';
end

Res.prm.N = size(M,1);
Res.prm.Ns = size(M,2);

if any(MinSample)
    Res.prm.MinSample = MinSample;
else 
    Res.prm.MinSample = 100;
end
if any(SampSizePerc)
    Res.prm.SampleSizesPerc = SampSizePerc;
else
    Res.prm.SampleSizesPerc = [.7 .7875 .9 1];
end
if any(Trials)
    Res.prm.Trials = Trials;
else
    Res.prm.Trials = [21 16 12 1];
end
if any(b_star)
    Res.prm.b_star = b_star;
else
    Res.prm.b_star = floor(.5*sqrt(Res.prm.Ns));  % Notice - this default value assumes no missing-values in M.
end
if any(MISS_VAL)    
    if any(find(ismember(1:Res.prm.b_star,MISS_VAL)))
        error (sprintf('Please use a MISS_VAL value which is not among 1,2,...,%d',Res.prm.b_star));
    end
    Res.prm.MISS_VAL = MISS_VAL;
else
    Res.prm.MISS_VAL = inf;
end

% Default parameters
Res.prm.RunSeed = 0;
rand ('state',Res.prm.RunSeed);
Res.prm.RunStart = datestr(now);
Res.prm.Documentation = 'Produced by Find_MI_relations Version 1.0';

Res.I = zeros(Res.prm.N);
Res.SampSize = zeros(Res.prm.N);

for n=1:Res.prm.N
    Res.GoodInds{n} = find(M(n,:)~=Res.prm.MISS_VAL);
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = QuantizeM (M,Res)

for binnum = 2 : Res.prm.b_star
    fprintf ('Quantize in advance each of the %d rows into %d bins ...\n',Res.prm.N,binnum);
    tmpQM = zeros(size(M));
    for n = 1 : Res.prm.N
        good_inds = Res.GoodInds{n}; 
        v = M(n,good_inds);
        [qv,cutoff] = MaxEntQuantize_v (v,binnum);
        v = M(n,:);
        v(good_inds) = qv;
        tmpQM(n,:) = v; % A quantized version of v into binnum bins where the MISSING VALS INCLUDED
    end
    Res.QM{binnum} = tmpQM;
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = FindPairInfo (ind1,ind2,Res)        

b_star = Res.prm.b_star;
SampleSizesPerc = Res.prm.SampleSizesPerc;
Trials = Res.prm.Trials;
MISS_VAL = Res.prm.MISS_VAL;

EstI = ones(1,b_star)*(-inf);
EstI_std = zeros(1,b_star);
BestEstI = 0;
BestEstI_binnum = 0;

% Find joint indices which are not missing-values
qv1 = Res.QM{2}(ind1,:);    
good_inds1 = Res.GoodInds{ind1};    
qv2 = Res.QM{2}(ind2,:);        
qv2 = qv2(good_inds1);    
good_inds2 = find(qv2 ~= MISS_VAL);    
good_inds = good_inds1(good_inds2);
sampsize = length(good_inds);
Res.SampSize(ind1,ind2) = sampsize;
Res.SampSize(ind2,ind1) = sampsize;

if sampsize < Res.prm.MinSample % Too small sample size
    Res.I(ind1,ind2) = -inf;
    Res.I(ind2,ind1) = -inf;
    return
end

for binnum = 2 : b_star

    qv1 = Res.QM{binnum}(ind1,:);
    qv2 = Res.QM{binnum}(ind2,:);

    % Leave only elements for which both v1 and v2 do NOT have a missing value
    qv1 = qv1(good_inds);
    qv2 = qv2(good_inds);
    
    if Res.prm.RandMode
        perm = randperm(qv_length);
        qv2 = qv2(perm);
    end
    
    [Info{binnum}] = DirectFindInfo (qv1,qv2,binnum,binnum,SampleSizesPerc,Trials);

    EstI (binnum) = Info{binnum}.EstI;
    EstI_std (binnum) = Info{binnum}.EstI_std;
    min_tmp_EstI = EstI(binnum) - .5*EstI_std(binnum); % reduce the error-bar from the current estimate
    max_prev_EstI = max ( EstI(1:binnum-1) + .5*EstI_std(1:binnum-1) ); % max of previous values + their errorbars
    if min_tmp_EstI > max_prev_EstI % improvement beyond the errorbars
        BestEstI = EstI (binnum);
        BestEstI_binnum = binnum;
    end
    
end

Res.I(ind1,ind2) = BestEstI;
Res.I(ind2,ind1) = BestEstI;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Quantize v into binnum bins which are -- more or less -- equally
% populated.
function [qv,cutoff] = MaxEntQuantize_v (v,binnum)

if isempty(v) % only missing values 
    qv = [];
    cutoff = [];
    return
end

qv = v;
[sv sv_inds] = sort(v);
group_size1 = floor(length(v)/binnum);
group_size2 = ceil(length(v)/binnum);   
extra = length(v) - binnum*group_size1;
cutoff = zeros(1,binnum-1);

% Randomly choose extra bins that will get another one point
bin_group_size = rand(1,binnum);
[tmp sinds] = sort(bin_group_size);
bin_group_size(sinds(1:(binnum-extra)))=1;
bin_group_size(sinds(end-extra+1:end))=0;

cutoff_ind = 0;
for b=1:binnum-1
    if bin_group_size(b)
        cutoff_ind = cutoff_ind + group_size1;
    else
        cutoff_ind = cutoff_ind + group_size2;
    end
    cutoff(b) = sv(cutoff_ind);
    if b==1
        b_inds = find( v <= cutoff(b) );
    else
        b_inds = find( (v > cutoff(b-1)) & (v <= cutoff(b)) );
    end
    qv(b_inds) = b;
end
b_inds = find(v > cutoff(binnum-1));
qv(b_inds) = binnum;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This generic function is used to estimate the Mutual Information (MI) between
% two already QUANTIZED vectors of the same dimension, via the direct
% estimation method. 
function [InfoOut] = DirectFindInfo (qv1,qv2,D1,D2,SampleSizesPerc,Trials)

FullSample = length(qv1);
SampleSizes = ceil(SampleSizesPerc*FullSample);
N = length(SampleSizes);
Ivals = zeros(N,max(Trials));
meanIvals = zeros(N,1);

for n=1:N % Loop over all the different Sample Sizes     
    tmpSampleSize = SampleSizes(n);    
    tmpv_leng = tmpSampleSize;            
    for t=1:Trials(n) % Loop over all the different trials for a specific sample size    
        if tmpSampleSize < FullSample % Choose a sample of v1 and v2 values
            perm = randperm(FullSample);             
            inds = perm(1:tmpSampleSize);            
            tmpqv1 = qv1(inds);            
            tmpqv2 = qv2(inds);            
        else
            tmpqv1 = qv1; 
            tmpqv2 = qv2;
        end            
        tmpCounts = ConstructCountMatrix (tmpqv1,tmpqv2,D1,D2); % Construct a count matrix out of the quantized vectors                
        Ivals(n,t) = FindInfoLocalMI (tmpCounts); % Find EMPIRICAL MI for the sample taken from v1 and v2        
    end % over Trials
    meanIvals(n) = mean(Ivals(n,1:Trials(n)));
end % over different Sample Sizes

[FitModel,Goodness,Output] = fit (1./SampleSizes',meanIvals,'poly1'); % linear extrapolation
Extrap_I = FitModel.p2; % value of the linear fit when x=0 (i.e., where the sample size goes to \infty).
    
InfoOut.Extrap_I = Extrap_I;
InfoOut.EstI = mean(Extrap_I);
[tmp MinSampSizeInd] = min(SampleSizes); % This should usually be equal to 1, but just in case.
InfoOut.EstI_std = std( Ivals(MinSampSizeInd,:) ); % take the std of the I vals obtained for the minimal sample size.

InfoOut.FitModel = FitModel;
InfoOut.Goodness = Goodness;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Counts = ConstructCountMatrix (qv1,qv2,D1,D2)

Counts = zeros(D1,D2);        
for l1=1:D1
    inds_l1 = find(qv1==l1);
    for l2=1:D2
        inds_l1l2 = find(qv2(inds_l1)==l2);
        Counts(l1,l2) = length(inds_l1l2);
    end 
end  

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MI = FindInfoLocalMI (tmpC)

tmpC = tmpC+eps;
sumC = sum(sum(tmpC));
Pxy = tmpC./sumC;
Px=sum(Pxy');
Hx=-sum(Px.*LocalLog(Px));
Py=sum(Pxy);
Hy=-sum(Py.*LocalLog(Py));

MI = Hx + Hy + sum(sum(Pxy.*LocalLog(Pxy)));

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = LocalLog (inp)

out = log2(inp);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
