function plotYpath( ypath)

global plotNumber T

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(ypath,'LineWidth',2)
hold on
hold on
title('Time path of income')


end

