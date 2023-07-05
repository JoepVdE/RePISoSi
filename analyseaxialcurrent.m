%sum(IArrayArchive(:,numPoints:end),2)
j = 1;
figure(2)
for i = numPoints:length(x1Array)


plotarrayx(j) = x1Array(i);
plotarrayx(j+1) = x2Array(i);
plotarrayy(j) = y1Array(i);
plotarrayy(j+1) = y2Array(i);
plotarrayz(j) = z1Array(i);
plotarrayz(j+1) = z2Array(i);
  j = j+2;

end

plot3(x1Array,y1Array,z1Array,'o','Color','b','MarkerSize',1)
hold on
plot3(plotarrayx,plotarrayy,plotarrayz,'-o','Color','b','MarkerSize',10,'MarkerFaceColor','r')
hold off