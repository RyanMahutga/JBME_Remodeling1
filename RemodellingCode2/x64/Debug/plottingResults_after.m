% errorlineting averages
close all
clear all
clc
net_num=1;
str2 = {'A0C_','A10C_','A25C_','A50C_','A100C_'};
str3 = {'pE_pC_','pE_npC_','npE_pC_'};
str4={'O','OS','U','S','US','NS'};
for m = 1:4
    len1=2;
    for n=2:len1
        for p=1:4
            w=p+1;
            
            str = strcat(str2{w},str3{2},num2str(net_num));
            
            old=cd(str);
            
            t = load(strcat('time',str4{m},'.txt'));
            time = t(2:end);
            
            orient(:,:,p) = load(strcat('orientation',str4{m},'.txt'));
            %pressurep(:,p) = load('pressure.txt');
            
            %fib_stress0 = load('fib_stress0.txt');
            %fib_stress = load('fib_stress5000.txt');
            %fib_type = load('fib_type5000.txt');
            %fib_rads = load('fib_rads5000.txt');
            %init_lens = load('init_lens5000.txt');
            
            %nodes = load('nodeData699.txt');
            %  figure
            % plot3(nodes(2:end,1),nodes(2:end,2), nodes(2:end,3),'ko')
            
            phi(:,:,p)=load(strcat('phi',str4{m},'.txt'));
            
            tot_phi(:,p) = phi(:,1,p) + phi(:,2,p) + phi(:,3,p);
            
            stretch(:,:,p)= load(strcat('stretch',str4{m},'.txt'));
            
            vol(:,p) = stretch(:,1,p).*stretch(:,2,p).*stretch(:,3,p);
            
            xstress(:,:,p) = load(strcat('xstress',str4{m},'.txt'))./1000;
            ystress(:,:,p) = load(strcat('ystress',str4{m},'.txt'))./1000;
            zstress(:,:,p) = load(strcat('zstress',str4{m},'.txt'))./1000;
            xzstress(:,:,p) = load(strcat('xzstress',str4{m},'.txt'))./1000;
            
            %         if p==1
            %            orient(:,:,n) = mean(orientp,3);
            % %            std_orient(:,:,n) = std(orientp,0,3);
            %            std_orient(:,:,n) = max(tinv([0.025 0.975],p-1))*std(orientp,0,3)/sqrt(p);
            %
            %            phi(:,:,n) = mean(phip,3);
            % %            std_phi(:,:,n) = std(phip,0,3);
            %            std_phi(:,:,n) = max(tinv([0.025 0.975],p-1))*std(phip,0,3)/sqrt(p);
            %
            %            stretch(:,:,n) = mean(stretchp,3);
            % %            std_stretch(:,:,n) = std(stretchp,0,3);
            %            std_stretch(:,:,n) = max(tinv([0.025 0.975],p-1))*std(stretchp,0,3)/sqrt(p);
            %
            %            vol(:,n) = mean(volp,2);
            % %            std_vol(:,n) = std(volp,0,2);
            %            std_vol(:,n) = max(tinv([0.025 0.975],p-1))*std(volp,0,2)/sqrt(p);
            %
            %            xstress(:,:,n) = mean(xstressp,3);
            % %            std_xstress(:,:,n) = std(xstressp,0,3);
            %            std_xstress(:,:,n)= max(tinv([0.025 0.975],p-1))*std(xstressp,0,3)/sqrt(p);
            %
            %            ystress(:,:,n) = mean(ystressp,3);
            % %            std_ystress(:,:,n) = std(ystressp,0,3);
            %            std_ystress(:,:,n)= max(tinv([0.025 0.975],p-1))*std(ystressp,0,3)/sqrt(p);
            %
            %            zstress(:,:,n) = mean(zstressp,3);
            % %            std_zstress(:,:,n) = std(zstressp,0,3);
            %            std_zstress(:,:,n)= max(tinv([0.025 0.975],p-1))*std(zstressp,0,3)/sqrt(p);
            %
            %            xystress(:,:,n) = mean(xystressp,3);
            % %            std_xystress(:,:,n) = std(xystressp,0,3);
            %            std_xystress(:,:,n)= max(tinv([0.025 0.975],p-1))*std(xystressp,0,3)/sqrt(p);
            %         end
            %
            %         fib_data0 = load('fiberData0.txt');
            %         fib_dataE = load('fiberData499.txt');
            %
            %         fibers0 = fib_data0(:,1:5);
            %         fibers0(:,1)=fibers0(:,1)+1;
            %         fibers0(:,2)=fibers0(:,2)+1;
            %
            %         fibers = fib_dataE(:,1:5);
            %         fibers(:,1)=fibers(:,1)+1;
            %         fibers(:,2)=fibers(:,2)+1;
            %
            %         fib_type = fib_data0(:,6);
            %
            %         net_data0 = load('nodeData0.txt');
            %         net_data = load('nodeData499.txt');
            %
            %         num_nodes = net_data0(1,1);
            %
            %         num_bnd0 = net_data0(1,2);
            %         num_bnd = net_data(1,2);
            %
            %         nodes0 = net_data0(2:num_nodes+1,:);
            %         nodes = net_data(2:num_nodes+1,:);
            %
            %         bnd_nodes0 = net_data0(num_nodes+2:1+num_nodes+num_bnd0,:);
            %         bnd_nodes = net_data(num_nodes+2:1+num_nodes+num_bnd,:);
            %
            %         fib_rads = fib_dataE(:,7);
            %         init_lens = fib_dataE(:,9);
            %         fib_stress = fib_dataE(:,12);
            %         fib_stretch = fib_dataE(:,10);
            %
            %         fib_rads0 = fib_data0(:,7);
            %         init_lens0 = fib_data0(:,9);
            %         fib_stress0 = fib_data0(:,12);
            %         fib_stretch0 = fib_data0(:,10);
            %
            %         crosscol = fibers(:,3) ;
            %         noncross = fib_rads(fib_type ==3 & crosscol == 0);
            %         cross = (fib_rads(fib_type ==3 & crosscol ~= 0));
            %
            %         histogram(noncross)
            %         hold on
            %         if n==2
            %             histogram(cross)
            %         end
            %
            %         lens0 = init_lens0(fib_type==1);
            %         rads0=fib_rads0(fib_type==3);
            %         rads=fib_rads(fib_type==3);
            %         lens = init_lens(fib_type==1);
            %         rads01=fib_rads0(fib_type==1);
            %         rads1=fib_rads(fib_type==1);
            %         lens03 = init_lens0(fib_type==3);
            %         lens3 = init_lens(fib_type==3);
            %
            %         fibers10 = fibers0(fib_type==3,:);
            %         fibers1 = fibers(fib_type==3,:);
            %
            cd(old)
        end
        % R0(:,:,n) = calc_periodic_orient(nodes, fibers10, rads01,lens0);
        % R(:,:,n) = calc_periodic_orient(nodes, fibers1,rads1,lens);
        %
        % val1 = eig(R0(:,:,n));
        % val2 = eig(R(:,:,n));
        % align(1) = max(val1)/min(val1);
        % align(2) = max(val2)/min(val2);
        
        % figure
        % t=pie([phi(1,1,n),phi(1,2,n),phi(1,3,n),1-sum(phi(1,:,n))])
        % colormap([1 0.7 0.2; 0 0 0; 0.7 0 0; 1 1 1])
        % t(2).FontSize = 16;
        % t(4).FontSize = 16;
        % t(6).FontSize = 16;
        % t(8).FontSize = 16;
        %
        % figure
        % t=pie([phi(end,1,n),phi(end,2,n),phi(end,3,n),1-sum(phi(end,:,n))])
        % colormap([1 0.7 0.2; 0 0 0; 0.7 0 0; 1 1 1])
        % t(2).FontSize = 16;
        % t(4).FontSize = 16;
        % t(6).FontSize = 16;
        % t(8).FontSize = 16;
        
    end
    
    % blue is unconnected -> red is fully connected
    xmin=-5;
    xmax=115;
    len2=4;
    % compare actin connectivity
    figure
    subplot(2,2,1)
    hold on
    for n=1:len2
        c = [(n-1)/5, 0, 1/n];
        plot(time,xstress(:,1,n),'LineStyle','-','Color', c,'LineWidth',2)
        plot(time,xstress(:,2,n),'LineStyle','--','Color', c,'LineWidth',2)
        plot(time,xstress(:,3,n),'LineStyle','-.','Color', c,'LineWidth',2)
        plot(time,xstress(:,4,n),'LineStyle',':','Color', c,'LineWidth',2)
    end
    legend('Total','Actin','Elastin','Collagen')
    xlabel('Time [h]')
    ylabel('Stress [kPa]')
    set(gca,'FontSize',16)
    axis([xmin xmax -10 350])
    
%     figure
subplot(2,2,2)
    hold on
    for n=1:len2
        c = [(n-1)/5, 0, 1/n];
        plot(time,xzstress(:,1,n),'LineStyle','-','Color', c,'LineWidth',2)
        plot(time,xzstress(:,2,n),'LineStyle','--','Color', c,'LineWidth',2)
        plot(time,xzstress(:,3,n),'LineStyle','-.','Color', c,'LineWidth',2)
        plot(time,xzstress(:,4,n),'LineStyle',':','Color', c,'LineWidth',2)
    end
    legend('Total','Actin','Elastin','Collagen')
    xlabel('Time [h]')
    ylabel('Shear Stress [kPa]')
    set(gca,'FontSize',16)
    axis([xmin xmax -5 45])
    
%     figure
subplot(2,2,3)
    hold on
    for n=1:len2
        c = [(n-1)/5, 0, 1/n];
        plot(time,vol(:,n),'LineStyle','-','Color', c,'LineWidth',2)
        plot(time,stretch(:,1,n),'LineStyle','--','Color', c,'LineWidth',2)
        plot(time,stretch(:,2,n),'LineStyle','-.','Color', c,'LineWidth',2)
        plot(time,stretch(:,3,n),'LineStyle',':','Color', c,'LineWidth',2)
    end
    legend('Total Vol.','Circ.','Axial','Radial')
    xlabel('Time [h]')
    ylabel('Stretch or Volume')
    set(gca,'FontSize',16)
    axis([xmin xmax 0.5 3.5])
    
%     figure
subplot(2,2,4)
    hold on
    for n=1:len2
        c = [(n-1)/5, 0, 1/n];
        plot(time,tot_phi(:,n),'LineStyle','-','Color', c,'LineWidth',2)
        plot(time,phi(:,1,n),'LineStyle','--','Color', c,'LineWidth',2)
        plot(time,phi(:,2,n),'LineStyle','-.','Color', c,'LineWidth',2)
        plot(time,phi(:,3,n),'LineStyle',':','Color', c,'LineWidth',2)
    end
    legend('Total','Actin','Elastin','Collagen')
    xlabel('Time [h]')
    ylabel('Fiber Volume Fraction, \phi')
    set(gca,'FontSize',16)
    axis([xmin xmax 0 1.0])
    
%     figure
%     hold on
%     for n=1:len1
%         yyaxis left
%         c = [(n-1)/5, 0, 1/n];
%         plot(time,orient(:,1,n),'LineStyle','--','Color', c,'LineWidth',2)
%         plot(time,orient(:,2,n),'LineStyle','-.','Color', c,'LineWidth',2)
%         plot(time,orient(:,3,n),'LineStyle',':','Color', c,'LineWidth',2)
%         yyaxis right
%         plot(time,orient(:,4,n),'LineStyle','-','Color', c,'LineWidth',2)
%     end
%     
%     yyaxis left
%     ylabel('Orientation')
%     ylim([0 1])
%     
%     yyaxis right
%     ylabel('Eigenvalue Ratio, \lambda_m_a_x / \lambda_m_i_n')
%     ylim([10 60])
%     
%     xlabel('Time [h]')
%     set(gca,'FontSize',16)
%     xlim([xmin xmax])
end

