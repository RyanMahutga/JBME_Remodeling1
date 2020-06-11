
time = [1 18 60 240 480 720];
Col = [5.9 14.4 17.2 33.9 31.8 36.6];
El = [7.7 27.1 23.3 18.2 15.9 9.5];
Cell = 0.45.*[44.9 22.2 29.3 24.9 24.7 25.4]; % cell is 45% fiber

ColS = [0.2 4.4 4.1 4.9 5.4 1.9];
ElS=[3.6 8.1 4.2 4.3 1.7 3.7];
CellS = 0.45.*[11.5 4.2 4.2 0.9 4.4 4.1];

figure
errorbar(time,Cell,CellS ,'Marker','o','LineStyle','--','Color',[0.9 0.7 0.1],'LineWidth',2)
hold on
errorbar(time,El,ElS ,'Marker','o','LineStyle','--','Color',[0 0 0],'LineWidth',2)
errorbar(time,Col,ColS ,'Marker','o','LineStyle','--','Color',[0.7 0.1 0.2],'LineWidth',2)
errorbar(time, Cell+Col+El, CellS+ElS+ColS,'Marker','o','LineStyle','--','Color',[0 0 0.9],'LineWidth',2)

xlabel('Time [days]')
ylabel('Fiber Volume Fraction, \phi [%]')
set(gca,'FontSize',16)
axis([-30 750 0 80])