function attach= dettach(r,attach)

attach(r)=0; % dettach cell r from whatever it was attached to
stop = 0;
dettach_ind = find(attach==r);   % find what is attached to cell r
if isempty(dettach_ind)
    stop=1;                     % if nothing is attached, no nothing more
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
end