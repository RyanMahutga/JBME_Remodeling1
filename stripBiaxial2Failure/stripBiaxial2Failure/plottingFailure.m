% boundedlineting averages
close all
clear all
clc
net_num=1;

fstr = "_sf";


failed = load(strcat('failed',fstr,'.txt'));

failed(failed < 0) = NaN;

phi=load(strcat('phi',fstr,'.txt'));

tot_phi = phi(:,1) + phi(:,2) + phi(:,3);

stretch= load(strcat('stretch',fstr,'.txt'));
stretch0 = stretch(1,:);
start_stretch = stretch(3,3);

stretch(:,1) = stretch(:,1)./stretch(1,1);
stretch(:,2) = stretch(:,2)./stretch(1,2);
stretch(:,3) = stretch(:,3)./stretch(1,3);

vol = stretch(:,1).*stretch(:,2).*stretch(:,3);

xstress = load(strcat('xstress',fstr,'.txt'))./1000;
ystress= load(strcat('ystress',fstr,'.txt'))./1000;
zstress = load(strcat('zstress',fstr,'.txt'))./1000;

xzstress = load(strcat('xzstress',fstr,'.txt'))./1000;
xzstress(xzstress<-10)=NaN;
xzstress(xzstress>10000)=NaN;
xstress(xstress<-10)=NaN;

if fstr == "_sf"
    shear =  load('shear_sf.txt');
    [~,idx_sh] = max(shear);
    shear(shear < 0) = NaN;
    
    angles = atand(shear./start_stretch);
    
    [pks_shear,locs_shear,widths_shear,proms_shear] = findpeaks(xzstress(:,1));
    
    idxpeaks = find(widths_shear>1.5,1);
    
    y_stretch_s = angles(locs_shear(idxpeaks));
    y_stress_s = xzstress(locs_shear(idxpeaks),1);
    
    x = angles;
    y=xzstress;
    x=x(x>0.05 & x<0.95);
    y=y(x>0.05 & x<0.95);
    
else
    
    [~,idx_s] = max(stretch(:,1) < 0);
    stretch(stretch < 0) = NaN;
    
    [pks_strip,locs_strip,widths_strip,proms_strip] = findpeaks(xstress(1:idx_s,3));
    
    idxpeaks = find(proms_strip>20,1);
    
    y_stretch_s = stretch(locs_strip(idxpeaks),1);
    y_stress_s = xstress(locs_strip(idxpeaks),1);
    
end

if fstr == "_sf"
    time = angles;
    xmin= -0.1;
    xmax = 3.0;
    stretchstr = ['Shear Angle [',char(176),']'] ;
    stressstr= 'Shear Stress [kPa]';
else
    time = stretch(:,1);
    xmin = 0.95;
    xmax = 2.25;
    stretchstr = 'Circ. Stretch, \lambda_1_1';
    stressstr= 'Circ. Stress [kPa]';
end

if fstr == "_f"
    figure
    plot(stretch (:,1),xstress(:,1),'o-','Color',[0 0 0.8], 'LineWidth',2)
    hold on
    plot(stretch (:,1),xstress(:,2),'o-','Color',[0.8 0.75 0.05], 'LineWidth',2)
    plot(stretch (:,1),xstress(:,3),'o-','Color',[0 0 0], 'LineWidth',2)
    plot(stretch (:,1),xstress(:,4),'o-','Color',[0.8 0 0], 'LineWidth',2)
    
    xlabel(stretchstr)
    ylabel(stressstr)
    legend('Total','Actin','Elastin','Collagen')
    
    set(gca,'FontSize',16)
    
    figure
    plot(stretch (:,1),ystress(:,1),'o-','Color',[0 0 0.8], 'LineWidth',2)
    hold on
    plot(stretch (:,1),ystress(:,2),'o-','Color',[0.8 0.75 0.05], 'LineWidth',2)
    plot(stretch (:,1),ystress(:,3),'o-','Color',[0 0 0], 'LineWidth',2)
    plot(stretch (:,1),ystress(:,4),'o-','Color',[0.8 0 0], 'LineWidth',2)
    
    xlabel(stretchstr)
    ylabel('Axial Stress, [kPa]')
    legend('Total','Actin','Elastin','Collagen')
    
    set(gca,'FontSize',16)
    
    figure
    plot(stretch (:,1),failed(:,1),'o-','Color',[0.8 0.75 0.05], 'LineWidth',2)
    hold on
    plot(stretch (:,1),failed(:,2),'o-','Color',[0 0 0], 'LineWidth',2)
    plot(stretch (:,1),failed(:,3),'o-','Color',[0.8 0 0], 'LineWidth',2)
    
    xlabel(stretchstr)
    ylabel('# of Failed Fibers')
    legend('Actin','Elastin','Collagen')
    
    set(gca,'FontSize',16)
end

if fstr == "_sf"
    figure
    plot(shear(:,1),xzstress(:,1),'o-','Color',[0 0 0.8], 'LineWidth',2)
    hold on
    plot(shear(:,1),xzstress(:,2),'o-','Color',[0.8 0.75 0.05], 'LineWidth',2)
    plot(shear(:,1),xzstress(:,3),'o-','Color',[0 0 0], 'LineWidth',2)
    plot(shear(:,1),xzstress(:,4),'o-','Color',[0.8 0 0], 'LineWidth',2)
    
    xlabel(stretchstr)
    ylabel('Shear Stress, [kPa]')
    legend('Total','Actin','Elastin','Collagen')
    
    set(gca,'FontSize',16)
    
     figure
    plot(shear(:,1),xstress(:,1),'o-','Color',[0 0 0.8], 'LineWidth',2)
    hold on
    plot(shear(:,1),xstress(:,2),'o-','Color',[0.8 0.75 0.05], 'LineWidth',2)
    plot(shear(:,1),xstress(:,3),'o-','Color',[0 0 0], 'LineWidth',2)
    plot(shear(:,1),xstress(:,4),'o-','Color',[0.8 0 0], 'LineWidth',2)
    
    xlabel(stretchstr)
    ylabel('Circ. Stress, [kPa]')
    legend('Total','Actin','Elastin','Collagen')
    
    set(gca,'FontSize',16)
   
    
    figure
    plot(shear(:,1),failed(:,1),'o-','Color',[0.8 0.75 0.05], 'LineWidth',2)
    hold on
    plot(shear(:,1),failed(:,2),'o-','Color',[0 0 0], 'LineWidth',2)
    plot(shear(:,1),failed(:,3),'o-','Color',[0.8 0 0], 'LineWidth',2)
    
    xlabel(stretchstr)
    ylabel('# of Failed Fibers')
    legend('Actin','Elastin','Collagen')
    
    set(gca,'FontSize',16)
end
    
