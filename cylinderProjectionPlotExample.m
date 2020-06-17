y1=sind(lat);
figure; scatter(lon(:),y1(:),16,p1(:),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360])
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7])