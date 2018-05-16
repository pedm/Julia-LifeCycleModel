function plotPath( path, title )

global plotnumber 


titleforgraph = title;
plotnumber = plotnumber + 1;
figure(plotnumber)
plot(path)
title(titleforgraph)


end

