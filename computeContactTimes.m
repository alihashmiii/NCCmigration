function contactTimes = computeContactTimes(attachments,dT)
% computes the contact times in the NC cell migration model from a cell
% array of attachment vectors at each time point
% each attachment vector has a list of cell indices that each cell is
% attached to (=0 for no attachment)

if nargin<2
    dT = 1;
end

numCells = length(attachments{end});
numTimePoints = length(attachments);
contactTimes = [];

timeCtr = numTimePoints; % we count backwards in time
for cellCtr=1:numCells
    % if this is a new cell (final timepoint), start a new counter
    thisContactTime = 1;
    for timeCtr = numTimePoints:-1:2 % step backwards in time
        prevNumCells = length(attachments{timeCtr-1});
        if cellCtr<=prevNumCells
            if attachments{timeCtr}(cellCtr)~=0
                % check if cell attachment is same in previous time
                if attachments{timeCtr}(cellCtr)==attachments{timeCtr-1}(cellCtr)
                    thisContactTime = thisContactTime + 1; % increase contact time
                else
                    % save contact time
                    contactTimes = [contactTimes; thisContactTime];
                    % start a new contactTime
                    thisContactTime = 1;
                end
            end
        else % this cell didn't exist in the previous time-step
            % check if contactTime for first time-step of this cell needs to be saved
            if attachments{timeCtr}(cellCtr)~=0
                % save contact time
                contactTimes = [contactTimes; thisContactTime];
            end
            break % break out of time-loop to go to next cell
        end
    end
    % check if contactTime for first time-step needs to be saved
    if cellCtr<=length(attachments{1})&&attachments{1}(cellCtr)~=0
        % save contact time
        contactTimes = [contactTimes; thisContactTime];
    end
end
% check that the overall sums are correct
assert(sum(contactTimes)==nnz(vertcat(attachments{:})),'sum of contact times does not match')
% adjust contact times for sampling interval
contactTimes = contactTimes*dT;
end

