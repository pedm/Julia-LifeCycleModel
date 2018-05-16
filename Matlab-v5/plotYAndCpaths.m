function plotYAndCpaths( ypath, cpath )

global plotNumber

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(ypath(:,1),'r','LineWidth',2)
hold on;
plot(cpath(:,1),'g','LineWidth',2)
hold on;
legend('Income','Consumption');
xlabel('Age');
title('Time path of income and consumption - individual 1')


plotNumber = plotNumber + 1;
figure(plotNumber)
plot(ypath(:,2),'r','LineWidth',2)
hold on;
plot(cpath(:,2),'g','LineWidth',2)
hold on;
legend('Income','Consumption');
xlabel('Age');
title('Time path of income and consumption - individual 2')

end

