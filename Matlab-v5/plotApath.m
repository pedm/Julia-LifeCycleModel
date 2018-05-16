function plotApath( apath, borrowingconstraints )

global plotNumber T

plotNumber = plotNumber + 1;
figure(plotNumber)
plot(apath,'LineWidth',2)
hold on
plot(borrowingconstraints, '--r','LineWidth', 2)
hold on
legend('1', '2', 'Borrowing con.');
xlabel('Age')
title('Time path of assets')


end

