% boundedlineting averages
close all
clear all
clc
net_num=1;

num_samples=6;

str2 = {'A100C_'};
str3 = {'pE_pC_','pE_npC_'};

subfldrs = {'Overload','Underload'};

plot_sub = 0;

fstr = "_f";

for z=1:1
    for w=1:5
        q=w;
        for p = 1:num_samples
            
            old2 = cd(strcat('PaperNets',num2str(p)));
            
            str = strcat(str2{q},str3{z},num2str(p));
            old=cd(str);
            
            if plot_sub ~= 0
                old3 = cd(subfldrs{plot_sub});
            end
            
            failedp(:,:,p) = load(strcat('failed',fstr,'.txt'));
            
            failedp(failedp < 0) = NaN;
            
            phip(:,:,p)=load(strcat('phi',fstr,'.txt'));
            
            tot_phip(:,p) = phip(:,1,p) + phip(:,2,p) + phip(:,3,p);
            
            stretchp(:,:,p)= load(strcat('stretch',fstr,'.txt'));
            start_stretch = stretchp(3,3,p);
            
            stretchp(:,1,p) = stretchp(:,1,p)./stretchp(1,1,p);
            stretchp(:,2,p) = stretchp(:,2,p)./stretchp(1,2,p);
            stretchp(:,3,p) = stretchp(:,3,p)./stretchp(1,3,p);
            
            volp(:,p) = stretchp(:,1,p).*stretchp(:,2,p).*stretchp(:,3,p);
            
            xstressp(:,:,p) = load(strcat('xstress',fstr,'.txt'))./1000;
            ystressp(:,:,p) = load(strcat('ystress',fstr,'.txt'))./1000;
            zstressp(:,:,p) = load(strcat('zstress',fstr,'.txt'))./1000;
            xzstressp(:,:,p) = load(strcat('xzstress',fstr,'.txt'))./1000;
            xzstressp(xzstressp<-10)=NaN;
            xzstressp(xzstressp>10000)=NaN;
            xstressp(xstressp<-10)=NaN;
            
            if fstr == "_sf"
                shearp(:,p) =  load('shear_sf.txt');
                [~,idx_sh] = max(shearp);
                shearp(shearp < 0) = NaN;
                
                angle(:,p) = atand(shearp(:,p)./start_stretch);
                
                [pks_shear,locs_shear,widths_shear,proms_shear] = findpeaks(xzstressp(:,1,p));
                
                y_stretch_s(p) = angle(locs_shear(1),p);
                y_stress_s(p) = xzstressp(locs_shear(1),1,p);
                
            else
                
                [~,idx_s] = max(stretchp(:,1,p) < 0);
                stretchp(stretchp < 0) = NaN;
                
                [pks_strip,locs_strip,widths_strip,proms_strip] = findpeaks(xstressp(1:idx_s,3,p));
                
                idx_loc = locs_strip(1);
                
                y_stretch_s(p) = stretchp(idx_loc,1,p);
                y_stress_s(p) = xstressp(idx_loc,1,p);
                
            end
            
            if plot_sub ~= 0
                cd(old3)
            end
            
            cd(old)
            cd(old2)
        end
        
        y_stretch(q) = mean(y_stretch_s);
        std_y_stretch(q) = max(tinv([0.025 0.975],p-1))*std(y_stretch_s)/sqrt(p);
        
        y_stress(q) = mean(y_stress_s);
        std_y_stress(q) = max(tinv([0.025 0.975],p-1))*std(y_stress_s)/sqrt(p);
        
        if fstr == "_sf"
            time = angle;
            xmin= -1;
            xmax = 91;
            stretchstr = ['Shear Angle [',char(176),']'] ;
            stressstr= 'Shear Stress [kPa]';
        else
            time = stretchp(:,1,:);
            xmin = 0.95;
            xmax = 2.25;
            stretchstr = 'Circ. Stretch, \lambda_1_1';
            stressstr= 'Circ. Stress [kPa]';
        end
        
        trans = 0.3;
        cmap = flipud(parula(6));
        
        % compare actin connectivity
                figure
                hold on
                for n=1:6
                    c = cmap(n,:);
                    if fstr == "_sf"
                        h1=plot(time(:,n),xzstressp(:,1,n),'LineStyle','-','LineWidth',2.5,'Color',c);
                        h2=plot(time(:,n),xzstressp(:,2,n),'LineStyle','--','LineWidth',2.5,'Color',c);
                        h3=plot(time(:,n),xzstressp(:,3,n),'LineStyle','-.','LineWidth',2.5,'Color',c);
                        h4=plot(time(:,n),xzstressp(:,4,n),'LineStyle',':','LineWidth',2.5,'Color',c);
                    else
                        h1=plot(time(:,n),xstressp(:,1,n),'LineStyle','-','LineWidth',2.5,'Color',c);
                        h2=plot(time(:,n),xstressp(:,2,n),'LineStyle','--','LineWidth',2.5,'Color',c);
                        h3=plot(time(:,n),xstressp(:,3,n),'LineStyle','-.','LineWidth',2.5,'Color',c);
                        h4=plot(time(:,n),xstressp(:,4,n),'LineStyle',':','LineWidth',2.5,'Color',c);
                    end
                end
                errorbar(y_stretch(q),y_stress(q),-std_y_stress(q),std_y_stress(q),-std_y_stretch(q),std_y_stretch(q),'k^','LineWidth',2, 'MarkerSize',6)
                %         errorbar(f_stretch(q),f_stress(q),-std_f_stress(q),std_f_stress(q),-std_f_stretch(q),std_f_stretch(q),'ko','LineWidth',2, 'MarkerSize',6)
                legend([h1, h2, h3, h4],'Total','Actin','Elastin','Collagen')
                xlabel(stretchstr)
                ylabel(stressstr)
                set(gca,'FontSize',16)
                xlim([xmin xmax])
        %         %-100 2500])
    end
end

cmap = flipud(parula(6));

figure
hold on
for q=1:6
    errorbar(y_stretch(q),y_stress(q),-std_y_stress(q),std_y_stress(q),...
        -std_y_stretch(q),std_y_stretch(q),'Marker','^','MarkerFaceColor',cmap(q,:),'Color',cmap(q,:),'LineWidth',2, 'MarkerSize',6)
end
%axis([1.4 2.2 300 900])
h = colorbar('Ticks',[0 0.25 0.5 0.75 1],'TickLabels',[0 10 25 50 100]);
set( h, 'YDir', 'reverse' );
xlabel(stretchstr)
ylabel(stressstr)
set(gca,'FontSize',16)

% y_m95(:,1) = y_stretch';
% y_m95(:,2) = std_y_stretch';
% y_m95(:,3) = y_stress';
% y_m95(:,4) = std_y_stress';
% 
% if plot_sub == 0
%     mkr = 'ko--';
% elseif plot_sub ==1
%     mkr='k>--';
% elseif plot_sub==2
%     mkr='k^--';
% elseif plot_sub==3
%     mkr='ks--';
% else
%     mkr='kv--';
% end

% ystr = y_stress; %./y_stress(1);
% xdata = [0,10,25,50,100];
% %figure
% % hold on
% % plot(xdata,ystr,mkr,'LineWidth',2, 'MarkerSize',6)
% % axis([-5 105 0.95 1.35])
% xlabel("Connectivity, [%]")
% ylabel("Normalized Failure Stress, [\sigma_n_% / \sigma_0_%]")
% set(gca,'FontSize',16)
% 
% datax = [datax;xdata];
% datay = [datay;ystr];
%legend('Initial','Shear','Overload','Overload + Shear','Underload')

