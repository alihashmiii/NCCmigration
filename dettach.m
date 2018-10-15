function attach= dettach(r,attach)

attach(r)=0;
stop = 0;
dettach_ind = find(attach==r);   % find what was attached to cell r
if isempty(dettach_ind)
    stop=1;                     % if nothing was attached, then detach
else
    attach(dettach_ind)=0;      % otherwise dettach all those cells
end

while stop==0
    next_dettach_ind = [];
    for ind=1:length(dettach_ind)
        next_dettach_ind = [next_dettach_ind; find(attach==dettach_ind(ind))];
    end
    if isempty(next_dettach_ind)
        stop = 1;
    else
        attach(next_dettach_ind)=0;
        dettach_ind = next_dettach_ind;
    end
end