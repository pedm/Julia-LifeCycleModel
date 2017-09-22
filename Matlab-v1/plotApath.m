function plotApath( apath )

global plotNumber 

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(apath,'LineWidth',2)
xlabel('Age')
ylabel('Assets')
hold on
title('Time path of assets')


end

