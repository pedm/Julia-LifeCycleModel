function plotVpath( vpath )

global plotNumber T

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(vpath,'LineWidth',2)
xlabel('Age');
ylabel('Value');
title('Time path of value')


end

