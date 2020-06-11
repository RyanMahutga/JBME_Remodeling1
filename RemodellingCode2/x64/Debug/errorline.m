function errorline(xdata, means, stdevs, color, style)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

plot(xdata,means,'Color',color,'LineStyle',style,'LineWidth',2)
hold on
plot(xdata, means+stdevs,'Color',color,'LineStyle',style,'LineWidth',0.5)
plot(xdata, means-stdevs,'Color',color,'LineStyle',style,'LineWidth',0.5)

end

