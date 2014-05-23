function handle = plotCircles(positions,radius,color,opacity)
% plots filled circles, e.g. cells
% Linus J. Schumacher 23.05.2014
%
% position should be in format 2-by-n
% opacity is an alhpa value

if size(positions,1) ~= 2
    positions = positions';
end

th = (0:0.1:2*pi)'; % for plotting the cell circles

% plot all the circles at once using Tony's trick for boradcasting
handle = fill(radius*cos(th(:,ones(size(positions(1,:))))) + positions(ones(size(th)),:),...
    radius*sin(th(:,ones(size(positions(2,:))))) + positions(2*ones(size(th)),:),color);

if opacity<1
    set(handle,'EdgeColor','none','FaceAlpha',opacity)
end

end

