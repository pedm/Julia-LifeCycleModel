function plotYCAndApaths( ypath, cpath, apath )

global plotNumber

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(ypath(:,1),'r','LineWidth',2)
hold on;
plot(cpath(:,1),'g','LineWidth',2)
hold on;
plot(apath(:,1),'b','LineWidth',2)
hold on;
legend('Income','Consumption','Assets');
xlabel('Age');
title('Time path of income and consumption - individual 1')


plotNumber = plotNumber + 1;
figure(plotNumber)
plot(ypath(:,2),'r','LineWidth',2)
hold on;
plot(cpath(:,2),'g','LineWidth',2)
hold on;
plot(apath(:,2),'b','LineWidth',2)
legend('Income','Consumption','Assets');
xlabel('Age');
title('Time path of income and consumption - individual 2')

end

