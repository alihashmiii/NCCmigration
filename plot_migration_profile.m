% plot migration profile
followerCounts = squeeze(cellDistributions(paramCtr,:,2,:));
leaderCounts = squeeze(cellDistributions(paramCtr,:,1,:));
lostCounts = squeeze(cellDistributions(paramCtr,:,3,:));
% % divide by total number of cells in run to get relative counts
% totalNumCells = sum(followerCounts + leaderCounts + lostCounts,2);
% followerCounts = followerCounts./totalNumCells(:,ones(1,size(followerCounts,2)));
% leaderCounts = leaderCounts./totalNumCells(:,ones(1,size(leaderCounts,2)));
% lostCounts = lostCounts./totalNumCells(:,ones(1,size(lostCounts,2)));

caProfile = mean(squeeze(caDistribution(paramCtr,:,:)));
caGrad = diff(caProfile)./diff(xlat_save')*filolength;
caCutOffIndcs = caGrad./sqrt(caProfile(1:end-1)) < sensingAccuracy;
[AX,H1,H2] = plotyy(xBins,mean(followerCounts),...
    xlat_save,caProfile);
hold(AX(1));
H3 = area(AX(1),xBins,[mean(lostCounts)', mean(followerCounts)', mean(leaderCounts)']);
hold(AX(2));
% mark the points where gradient in concentration is below sensing
% threshold
plot(AX(2),xlat_save(caCutOffIndcs),caProfile(caCutOffIndcs),'x','Color',[0 0.5 0]);
set(H2,'LineStyle','--','Color',[0 0.5 0]), set(H3,'LineStyle','-','EdgeColor','b')
set(H3(3),'FaceColor','y'), set(H3(2),'FaceColor','w'), set(H3(1),'FaceColor','r')
ylim(AX(1),[0, 12]), xlim(AX(1),[0 800]), xlim(AX(2),[0 800]), ylim(AX(2),[0 0.4])
set(AX(2),'YTick',[0 1])
set(AX(1),'YTick',[0 3 6 9 12])
grid on