function plotCpath( cpath )

global plotNumber T

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(cpath,'LineWidth',2)
xlabel('Age');
title('Time path of consumption')


end

