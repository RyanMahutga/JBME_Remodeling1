% boundedlineting averages
close all
clear all
clc
net_num=1;

num_samples=5;
str2 = {'A0C_','A10C_','A25C_','A50C_','A100C_'};
str3 = {'pE_pC_','pE_npC_','npE_pC_'};

subfldrs = {'Shear','Overload','Overload_Shear','Underload'};

plot_sub = 1;

for n=2:2
    for w=1:5
        q=w;
        for p =1:num_samples
            
            old2 = cd(strcat('PaperNets',num2str(p)));
            
            str = strcat(str2{q},str3{n},num2str(p));
            old=cd(str);
            
            if plot_sub ~= 0
               old3 = cd(subfldrs{plot_sub}); 
            end
            
            t = load('time.txt');
            time = t(2:end);
            
            orientp(:,:,p) = load('orientation.txt');
            %pressurep(:,p) = load('pressure.txt');
            
            %fib_stress0 = load('fib_stress0.txt');
            %fib_stress = load('fib_stress5000.txt');
            %fib_type = load('fib_type5000.txt');
            %fib_rads = load('fib_rads5000.txt');
            %init_lens = load('init_lens5000.txt');
            
            %nodes = load('nodeData699.txt');
            %  figure
            % plot3(nodes(2:end,1),nodes(2:end,2), nodes(2:end,3),'ko')
            
            phip(:,:,p)=load('phi.txt');
            
            tot_phip(:,p) = phip(:,1,p) + phip(:,2,p) + phip(:,3,p);
            
            stretchp(:,:,p)= load('stretch.txt');
            
            volp(:,p) = stretchp(:,1,p).*stretchp(:,2,p).*stretchp(:,3,p);
            
            xstressp(:,:,p) = load('xstress.txt')./1000;
            ystressp(:,:,p) = load('ystress.txt')./1000;
            zstressp(:,:,p) = load('zstress.txt')./1000;
            xzstressp(:,:,p) = load('xzstress.txt')./1000;
            
            if plot_sub ~= 0
                cd(old3)
            end
   
            cd(old)
            cd(old2)
        end
        orient(:,:,q) = mean(orientp,3);
        std_orient(:,:,q) = max(tinv([0.025 0.975],p-1))*std(orientp,0,3)/sqrt(p);
        
        phi(:,:,q) = mean(phip,3);
        std_phi(:,:,q) = max(tinv([0.025 0.975],p-1))*std(phip,0,3)/sqrt(p);
        
        tot_phi(:,:,q) = mean(tot_phip,3);
        std_tot_phi(:,:,q) = max(tinv([0.025 0.975],p-1))*std(tot_phip,0,3)/sqrt(p);
        
        stretch(:,:,q) = mean(stretchp,3);
        std_stretch(:,:,q) = max(tinv([0.025 0.975],p-1))*std(stretchp,0,3)/sqrt(p);
        
        vol(:,q) = mean(volp,2);
        std_vol(:,q) = max(tinv([0.025 0.975],p-1))*std(volp,0,2)/sqrt(p);
        
        xstress(:,:,q) = mean(xstressp,3);
        std_xstress(:,:,q)= max(tinv([0.025 0.975],p-1))*std(xstressp,0,3)/sqrt(p);
        
        ystress(:,:,q) = mean(ystressp,3);
        std_ystress(:,:,q)= max(tinv([0.025 0.975],p-1))*std(ystressp,0,3)/sqrt(p);
        
        zstress(:,:,q) = mean(zstressp,3);
        std_zstress(:,:,q)= max(tinv([0.025 0.975],p-1))*std(zstressp,0,3)/sqrt(p);
        
        xzstress(:,:,q) = mean(xzstressp,3);
        std_xzstress(:,:,q)= max(tinv([0.025 0.975],p-1))*std(xzstressp,0,3)/sqrt(p);
        
        %                     fib_data0 = load('fiberData0.txt');
        %                     fib_dataE = load('fiberData499.txt');
        %
        %                     fibers0 = fib_data0(:,1:5);
        %                     fibers0(:,1)=fibers0(:,1)+1;
        %                     fibers0(:,2)=fibers0(:,2)+1;
        %
        %                     fibers = fib_dataE(:,1:5);
        %                     fibers(:,1)=fibers(:,1)+1;
        %                     fibers(:,2)=fibers(:,2)+1;
        %
        %                     fib_type = fib_data0(:,6);
        %
        %                     net_data0 = load('nodeData0.txt');
        %                     net_data = load('nodeData499.txt');
        %
        %                     num_nodes = net_data0(1,1);
        %
        %                     num_bnd0 = net_data0(1,2);
        %                     num_bnd = net_data(1,2);
        %
        %                     nodes0 = net_data0(2:num_nodes+1,:);
        %                     nodes = net_data(2:num_nodes+1,:);
        %
        %                     bnd_nodes0 = net_data0(num_nodes+2:1+num_nodes+num_bnd0,:);
        %                     bnd_nodes = net_data(num_nodes+2:1+num_nodes+num_bnd,:);
        %
        %                     fib_rads = fib_dataE(:,7);
        %                     init_lens = fib_dataE(:,9);
        %                     fib_stress = fib_dataE(:,12);
        %                     fib_stretch = fib_dataE(:,10);
        %
        %                     fib_rads0 = fib_data0(:,7);
        %                     init_lens0 = fib_data0(:,9);
        %                     fib_stress0 = fib_data0(:,12);
        %                     fib_stretch0 = fib_data0(:,10);
        %
        %                     crosscol = fibers(:,3) ;
        %                     noncross = fib_rads(fib_type ==3 & crosscol == 0);
        %                     cross = (fib_rads(fib_type ==3 & crosscol ~= 0));
        %
        %                     histogram(noncross)
        %                     hold on
        %                     if n==2
        %                         histogram(cross)
        %                     end
        %
        %                     lens0 = init_lens0(fib_type==1);
        %                     rads0=fib_rads0(fib_type==3);
        %                     rads=fib_rads(fib_type==3);
        %                     lens = init_lens(fib_type==1);
        %                     rads01=fib_rads0(fib_type==1);
        %                     rads1=fib_rads(fib_type==1);
        %                     lens03 = init_lens0(fib_type==3);
        %                     lens3 = init_lens(fib_type==3);
        %
        %                     fibers10 = fibers0(fib_type==3,:);
        %                     fibers1 = fibers(fib_type==3,:);
        %
    end
end

% blue is unconnected -> red is fully connected
xmin=-5;
xmax=95;

trans = 0.5;
cmap = flipud(parula(5));

% compare actin connectivity
figure
subplot(2,2,1)
hold on
boundedline(time,xstress(:,1,1),std_xstress(:,1,1),...
     time,xstress(:,1,2),std_xstress(:,1,2),...
     time,xstress(:,1,3),std_xstress(:,1,3),...
     time,xstress(:,1,4),std_xstress(:,1,4),...
     time,xstress(:,1,5),std_xstress(:,1,5),...
     time,xstress(:,2,1),std_xstress(:,2,1),...
     time,xstress(:,2,2),std_xstress(:,2,2),...
     time,xstress(:,2,3),std_xstress(:,2,3),...
     time,xstress(:,2,4),std_xstress(:,2,4),...
     time,xstress(:,2,5),std_xstress(:,2,5),...
     time,xstress(:,3,1),std_xstress(:,3,1),...
     time,xstress(:,3,2),std_xstress(:,3,2),...
     time,xstress(:,3,3),std_xstress(:,3,3),...
     time,xstress(:,3,4),std_xstress(:,3,4),...
     time,xstress(:,3,5),std_xstress(:,3,5),...
     time,xstress(:,4,1),std_xstress(:,4,1),...
     time,xstress(:,4,2),std_xstress(:,4,2),...
     time,xstress(:,4,3),std_xstress(:,4,3),...
     time,xstress(:,4,4),std_xstress(:,4,4),...
     time,xstress(:,4,5),std_xstress(:,4,5),...
     'cmap', [cmap;cmap;cmap;cmap], ':', 'transparency',trans)

      for n=1:5
          c = cmap(n,:);
          h1=plot(time,xstress(:,1,n),'LineStyle','-','LineWidth',2.5,'Color',c);
          h2=plot(time,xstress(:,2,n),'LineStyle','--','LineWidth',2.5,'Color',c);
          h3=plot(time,xstress(:,3,n),'LineStyle','-.','LineWidth',2.5,'Color',c);
          h4=plot(time,xstress(:,4,n),'LineStyle',':','LineWidth',2.5,'Color',c);
      end
legend([h1, h2, h3, h4],'Total','Actin','Elastin','Collagen')
xlabel('Time [h]')
ylabel('Circ. Stress [kPa]')
set(gca,'FontSize',16)
xlim([xmin xmax]) ;% -10 350])


% for n=1:5
%     c = [(n-1)/4, 0, 1/n];
%     boundedline(time,xstress(:,1,n),std_xstress(:,1,n), 'cmap', c,'-', 'transparency',trans)
%     boundedline(time,xstress(:,2,n),std_xstress(:,2,n), 'cmap', c,'--', 'transparency',trans)
%     boundedline(time,xstress(:,3,n),std_xstress(:,3,n), 'cmap', c,'-.', 'transparency',trans)
%     boundedline(time,xstress(:,4,n),std_xstress(:,4,n), 'cmap', c,':', 'transparency',trans)
% end


subplot(2,2,2)
hold on
boundedline(time,xzstress(:,1,1),std_xzstress(:,1,1),...
     time,xzstress(:,1,2),std_xzstress(:,1,2),...
     time,xzstress(:,1,3),std_xzstress(:,1,3),...
     time,xzstress(:,1,4),std_xzstress(:,1,4),...
     time,xzstress(:,1,5),std_xzstress(:,1,5),...
     time,xzstress(:,2,1),std_xzstress(:,2,1),...
     time,xzstress(:,2,2),std_xzstress(:,2,2),...
     time,xzstress(:,2,3),std_xzstress(:,2,3),...
     time,xzstress(:,2,4),std_xzstress(:,2,4),...
     time,xzstress(:,2,5),std_xzstress(:,2,5),...
     time,xzstress(:,3,1),std_xzstress(:,3,1),...
     time,xzstress(:,3,2),std_xzstress(:,3,2),...
     time,xzstress(:,3,3),std_xzstress(:,3,3),...
     time,xzstress(:,3,4),std_xzstress(:,3,4),...
     time,xzstress(:,3,5),std_xzstress(:,3,5),...
     time,xzstress(:,4,1),std_xzstress(:,4,1),...
     time,xzstress(:,4,2),std_xzstress(:,4,2),...
     time,xzstress(:,4,3),std_xzstress(:,4,3),...
     time,xzstress(:,4,4),std_xzstress(:,4,4),...
     time,xzstress(:,4,5),std_xzstress(:,4,5),...
     'cmap', [cmap;cmap;cmap;cmap], ':', 'transparency',trans)

      for n=1:5
          c = cmap(n,:);
          h1=plot(time,xzstress(:,1,n),'LineStyle','-','LineWidth',2.5,'Color',c);
          h2=plot(time,xzstress(:,2,n),'LineStyle','--','LineWidth',2.5,'Color',c);
          h3=plot(time,xzstress(:,3,n),'LineStyle','-.','LineWidth',2.5,'Color',c);
          h4=plot(time,xzstress(:,4,n),'LineStyle',':','LineWidth',2.5,'Color',c);
      end

% for n=1:5
%     c = [(n-1)/4, 0, 1/n];
%     boundedline(time,xzstress(:,1,n),std_xzstress(:,1,n), 'cmap', c,'-', 'transparency',trans)
%     boundedline(time,xzstress(:,2,n),std_xzstress(:,2,n), 'cmap', c,'--', 'transparency',trans)
%     boundedline(time,xzstress(:,3,n),std_xzstress(:,3,n), 'cmap', c,'-.', 'transparency',trans)
%     boundedline(time,xzstress(:,4,n),std_xzstress(:,4,n), 'cmap', c,':', 'transparency',trans)
% end

legend([h1, h2, h3, h4],'Total','Actin','Elastin','Collagen')
xlabel('Time [h]')
ylabel('Shear Stress [kPa]')
set(gca,'FontSize',16)
xlim([xmin xmax])
% axis([xmin xmax -3 40])

% figure
subplot(2,2,3)
hold on
boundedline(time,vol(:,1),std_vol(:,1),...
     time,vol(:,2),std_vol(:,2),...
     time,vol(:,3),std_vol(:,3),...
     time,vol(:,4),std_vol(:,4),...
     time,vol(:,5),std_vol(:,5),...
     time,stretch(:,1,1),std_stretch(:,1,1),...
     time,stretch(:,1,2),std_stretch(:,1,2),...
     time,stretch(:,1,3),std_stretch(:,1,3),...
     time,stretch(:,1,4),std_stretch(:,1,4),...
     time,stretch(:,1,5),std_stretch(:,1,5),...
     time,stretch(:,2,1),std_stretch(:,2,1),...
     time,stretch(:,2,2),std_stretch(:,2,2),...
     time,stretch(:,2,3),std_stretch(:,2,3),...
     time,stretch(:,2,4),std_stretch(:,2,4),...
     time,stretch(:,2,5),std_stretch(:,2,5),...
     time,stretch(:,3,1),std_stretch(:,3,1),...
     time,stretch(:,3,2),std_stretch(:,2,2),...
     time,stretch(:,3,3),std_stretch(:,2,3),...
     time,stretch(:,3,4),std_stretch(:,2,4),...
     time,stretch(:,3,5),std_stretch(:,2,5),...
     'cmap', [cmap;cmap;cmap;cmap], ':', 'transparency',trans)

      for n=1:5
          c = cmap(n,:);
          h1=plot(time,vol(:,n),'LineStyle','-','LineWidth',2.5,'Color',c);
          h2=plot(time,stretch(:,1,n),'LineStyle','--','LineWidth',2.5,'Color',c);
          h3=plot(time,stretch(:,2,n),'LineStyle','-.','LineWidth',2.5,'Color',c);
          h4=plot(time,stretch(:,3,n),'LineStyle',':','LineWidth',2.5,'Color',c);
      end
%       
% for n=1:5
%     c = [(n-1)/4, 0, 1/n];
%     boundedline(time,vol(:,n),std_vol(:,n), 'cmap', c,'-', 'transparency',trans)
%     boundedline(time,stretch(:,1,n),std_stretch(:,1,n), 'cmap', c,'--', 'transparency',trans)
%     boundedline(time,stretch(:,2,n),std_stretch(:,2,n), 'cmap', c,'-.', 'transparency',trans)
%    boundedline(time,stretch(:,3,n),std_stretch(:,3,n), 'cmap', c,':', 'transparency',trans)
% end

legend([h1,h2,h3,h4],'Total Vol.','Circ.','Axial','Radial')
xlabel('Time [h]')
ylabel('Stretch or Volume Change')
set(gca,'FontSize',16)
xlim([xmin xmax])
% axis([xmin xmax 0.25 2.0])

% figure
subplot(2,2,4)
hold on
boundedline(time,tot_phi(:,1),std_tot_phi(:,1),...
     time,tot_phi(:,2),std_tot_phi(:,2),...
     time,tot_phi(:,3),std_tot_phi(:,3),...
     time,tot_phi(:,4),std_tot_phi(:,4),...
     time,tot_phi(:,5),std_tot_phi(:,5),...
     time,phi(:,1,1),std_phi(:,1,1),...
     time,phi(:,1,2),std_phi(:,1,2),...
     time,phi(:,1,3),std_phi(:,1,3),...
     time,phi(:,1,4),std_phi(:,1,4),...
     time,phi(:,1,5),std_phi(:,1,5),...
     time,phi(:,2,1),std_phi(:,2,1),...
     time,phi(:,2,2),std_phi(:,2,2),...
     time,phi(:,2,3),std_phi(:,2,3),...
     time,phi(:,2,4),std_phi(:,2,4),...
     time,phi(:,2,5),std_phi(:,2,5),...
     time,phi(:,3,1),std_phi(:,3,1),...
     time,phi(:,3,2),std_phi(:,2,2),...
     time,phi(:,3,3),std_phi(:,2,3),...
     time,phi(:,3,4),std_phi(:,2,4),...
     time,phi(:,3,5),std_phi(:,2,5),...
     'cmap', [cmap;cmap;cmap;cmap], ':', 'transparency',trans)

      for n=1:5
          c = cmap(n,:);
          h1=plot(time,tot_phi(:,n),'LineStyle','-','LineWidth',2.5,'Color',c);
          h2=plot(time,phi(:,1,n),'LineStyle','--','LineWidth',2.5,'Color',c);
          h3=plot(time,phi(:,2,n),'LineStyle','-.','LineWidth',2.5,'Color',c);
          h4=plot(time,phi(:,3,n),'LineStyle',':','LineWidth',2.5,'Color',c);
      end
      
% for n=1:5
%     c = [(n-1)/4, 0, 1/n];
%     boundedline(time,tot_phi(:,n),std_tot_phi(:,n), 'cmap', c,'-', 'transparency',trans)
%     boundedline(time,phi(:,1,n),std_phi(:,1,n), 'cmap', c,'--', 'transparency',trans)
%     boundedline(time,phi(:,2,n),std_phi(:,2,n), 'cmap', c,'-.', 'transparency',trans)
%     boundedline(time,phi(:,3,n),std_phi(:,3,n), 'cmap', c,':', 'transparency',trans)
% end

legend([h1,h2,h3,h4],'Total','Actin','Elastin','Collagen')
xlabel('Time [h]')
ylabel('Fiber Volume Fraction, \phi')
set(gca,'FontSize',16)
xlim([xmin xmax])
% axis([xmin xmax 0 0.7])

hp4 = get(subplot(2,2,4),'Position');
h=colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.02  hp4(2)+hp4(3)*2.1],...
    'Ticks',[0 0.25 0.5 0.75 1],'TickLabels',[0 10 25 50 100]);
set( h, 'YDir', 'reverse' );
h.Label.String = '% Actin Connectivity';

%
% % figure
% hold on
% yyaxis left
% for n=1:5
%    c = [(n-1)/5, 0, 1/n];
%    plot(time,orient(:,1,n),'LineStyle','--','Color', 'cmap', c,'LineWidth',2)
%    plot(time,orient(:,2,n),'LineStyle','-.','Color', 'cmap', c,'LineWidth',2)
%    plot(time,orient(:,3,n),'LineStyle',':','Color', 'cmap', c,'LineWidth',2)
% end
% ylabel('Orientation')
% ylim([0 1])
%
% yyaxis right
% hold on
% for n=1:5
%    %c = [(n-1)/5, 0, 1/n];
%    plot(time,orient(:,4,n),'k-','LineWidth',2)
% end
% ylabel('Eigenvalue Ratio, \lambda_m_a_x / \lambda_m_i_n')
% ylim([15 85])
%
% xlabel('Time [h]')
% set(gca,'FontSize',16)
% xlim([xmin xmax])
%
%
