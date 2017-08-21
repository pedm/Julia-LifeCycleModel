function plotCpath( cpath )

global plotNumber

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(cpath,'LineWidth',2)
xlabel('Age')
ylabel('Consumption')
title('Time path of consumption')

end

