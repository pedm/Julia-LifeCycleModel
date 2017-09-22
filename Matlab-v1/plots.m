%% Plots
global plotNumber

%%
% GRAPH THE POLICY FUNCTION    

 plotNumber = plotNumber + 1;
    figure(plotNumber);
    plot(Agrid(whichyear, plotNode1:plotNodeLast),policyC(whichyear, plotNode1:plotNodeLast),'g','LineWidth',2)
    hold on;
    xlabel('Asset');
    ylabel('Policy function (consumption function)');
    legend('Consumption');
    title('Policy function (consumption function)')
 
%%
 % GRAPH THE VALUE FUNCTION 
  plotNumber = plotNumber + 1;
    figure(plotNumber);
    plot(Agrid(whichyear, plotNode1:plotNodeLast),val(whichyear, plotNode1:plotNodeLast),'g','LineWidth',2)
    hold on;
    xlabel('Asset');
    ylabel('Value');
    legend('Value function');
    title('Value function')

 
   %%
 % COMPARE ANALYTICAL AND NUMERICAL CONSUMPTION 

plotNumber = plotNumber + 1;

figure(plotNumber);
plot(Agrid(whichyear, plotNode1:plotNodeLast),policyC(whichyear, plotNode1:plotNodeLast),'r','LineWidth',2)
hold on;
plot(Agrid(whichyear, plotNode1:plotNodeLast),policyC_analytic(whichyear, plotNode1:plotNodeLast),'g','LineWidth',2)
hold on;
xlabel('Asset');
ylabel('Policy (consumption)');
legend('Numerical','Analytical');
title('Comparing numerical and analytical consumption functions')   
