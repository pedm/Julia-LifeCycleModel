global plotNumber                                   % for graphing


%%
% GRAPH THE POLICY FUNCTION    

 plotNumber = plotNumber + 1;
    figure(plotNumber);
    plot(Agrid(whichyear, plotNode1:plotNodeLast),policyC(whichyear, plotNode1:plotNodeLast, numPointsY),'g','LineWidth',2)
    hold on;
    plot(Agrid(whichyear, plotNode1:plotNodeLast),policyC(whichyear, plotNode1:plotNodeLast, 1),'r','LineWidth',2)
    hold on;
    xlabel('Asset');
    ylabel('Policy function (consumption function)');
    legend('Higest income', 'Lowest income');
    title('Policy function (consumption function)')
 
%%
 % GRAPH THE VALUE FUNCTION 
  plotNumber = plotNumber + 1;
    figure(plotNumber);
    plot(Agrid(whichyear, plotNode1:plotNodeLast),val(whichyear, plotNode1:plotNodeLast, numPointsY),'g','LineWidth',2)
    hold on;
    plot(Agrid(whichyear, plotNode1:plotNodeLast),val(whichyear, plotNode1:plotNodeLast, 1),'r','LineWidth',2)
    hold on;
    xlabel('Asset');
    ylabel('Value');
    legend('Val. cond. on highest inc.','Val. cond. on lowest inc.');
    title('Value function')

 