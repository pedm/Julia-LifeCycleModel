function plotCZoomedOut( cpath )

global plotNumber T


lowerboundforyaxis = min(cpath)*0.9;
 upperboundforyaxis = max(cpath)*1.1;
 
 plotNumber = plotNumber + 1;
 figure(plotNumber)
 plot(cpath, 'LineWidth',2)
 axis([1 T lowerboundforyaxis upperboundforyaxis])
 xlabel('Age')
 ylabel('Consumption')
 title('Time path of consumption')

end
