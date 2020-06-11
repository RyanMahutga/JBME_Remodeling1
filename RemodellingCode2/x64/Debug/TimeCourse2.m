
h = 35e-6; %m
Ro = 0.5*1771e-6; %m

MAP = 100; %mmHg

Ri = Ro-h;
R = Ri+0.5*h;
MaxP = 13.33/100*MAP;

for t=1:492
    
    Time(t) = t-1-12; 
    P(t) = MaxP*(1-exp(-0.04*(Time(t)+12)));
    BP_mmHg(t) = P(t)*100/13.33;
    
    sig_C(t) = P(t)*Ri^2/(Ro^2-Ri^2) + P(t)*Ri^2*Ro^2/(R^2*(Ro^2-Ri^2)) ;
    sig_A(t) = P(t)*Ri^2/(Ro^2-Ri^2) ;
    sig_R(t) = P(t)*Ri^2/(Ro^2-Ri^2) - P(t)*Ri^2*Ro^2/(R^2*(Ro^2-Ri^2)) ;
end
% 
figure
plot(Time, sig_C, 'k-','LineWidth',2)
hold on
plot(Time, sig_A, 'k--','LineWidth',2)
plot(Time, sig_R, 'k:','LineWidth',2)
% Cell
plot(Time, 0.63.*sig_C, 'm-','LineWidth',2)
plot(Time, 0.63.*sig_A, 'm--','LineWidth',2)
plot(Time, 0.63.*sig_R, 'm:','LineWidth',2)
% ECM
plot(Time, 0.37.*sig_C, 'c-','LineWidth',2)
plot(Time, 0.37.*sig_A, 'c--','LineWidth',2)
plot(Time, 0.37.*sig_R, 'c:','LineWidth',2)
% 
% xlabel('Time (days)')
% ylabel('Boundary Stress (kPa)')
% legend('Circumferential','Axial', 'Radial')
% set(gca,'FontSize',16)
% axis([-10 485 -20 350])
% 
% % from Kassab
% BP=[0; 30.8558558558559;40.7657657657658;43.2432432432433;...
%     57.6576576576577;66.4414414414415;79.7297297297298;...
%     80.1801801801802;78.8288288288289];
% Time2 = [-12; 1.45624891718642; 2.97900965869716;6.05200753638254;...
%     8.92471197158698;13.8960444819820; 16.9360706860707;...
%     23.8775554400555;29.8920434857935];
% 
% figure
% plot(Time, BP_mmHg, 'k-','LineWidth',2)
% hold on
% plot(Time2, BP, 'k^','LineWidth',2)
% 
% xlabel('Time (days)')
% % ylabel('Mean Arterial Pressure (mmHg)')
% legend('Fit','Huang, Y. et al. (2005)')
% set(gca,'FontSize',16)
% axis([-10 50 -10 110])