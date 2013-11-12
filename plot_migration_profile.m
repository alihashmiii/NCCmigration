% plot migration profile
followerCounts = squeeze(cellDistributions(paramCtr,:,2,:));
leaderCounts = squeeze(cellDistributions(paramCtr,:,1,:));
lostCounts = squeeze(cellDistributions(paramCtr,:,3,:));
[AX,H1,H2] = plotyy(xBins,mean(followerCounts),...
    xBins,mean(squeeze(caDistribution(paramCtr,:,:)))/numRepeats);
hold(AX(1));
H3 = area(AX(1),xBins,[mean(lostCounts)', mean(followerCounts)', mean(leaderCounts)']);
set(H2,'LineStyle',':','Color',[0 0.5 0]), set(H3,'LineStyle','-','EdgeColor','b')
set(H3(3),'FaceColor','y'), set(H3(2),'FaceColor','w'), set(H3(1),'FaceColor','r')
ylim(AX(1),[0, 7]), xlim(AX(1),[0 735]), xlim(AX(2),[0 735])
set(AX(2),'YTick',[0 1])
set(AX(1),'YTick',[0 2 4 6])
