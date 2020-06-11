close all
clear all

time = load('time.txt');

orient = load('orientation.txt');
% pressure = load('Pressure.txt');

%fib_stress0 = load('fib_stress0.txt');
%fib_stress = load('fib_stress5000.txt');
%fib_type = load('fib_type5000.txt');
%fib_rads = load('fib_rads5000.txt');
%init_lens = load('init_lens5000.txt');

phi=load('phi.txt');
WSS=load('WSS.txt');

% fib_data0 = load('fiberData0.txt');
% % fib_dataE = load('fiberData254.txt');
% 
% fibers0 = fib_data0(:,1:5);
% fibers = fib_dataE(:,1:5);
% 
% fib_type = fib_data0(:,6);
% 
% net_data0 = load('nodeData0.txt');
% net_data = load('nodeData2.txt');
% 
% num_nodes = net_data0(1,1);
% 
% num_bnd0 = net_data0(1,2);
% num_bnd = net_data(1,2);
% 
% nodes0 = net_data(2:num_nodes+1,:);
% nodes = net_data(2:num_nodes+1,:);
% 
% bnd_nodes0 = net_data(num_nodes+2:1+num_nodes+num_bnd0,:);
% bnd_nodes = net_data(num_nodes+2:1+num_nodes+num_bnd,:);
% 
% fib_rads = fib_dataE(:,7);
% init_lens = fib_dataE(:,9);
% fib_stress = fib_dataE(:,12);
% fib_stretch = fib_dataE(:,10);
% 
% fib_rads0 = fib_data0(:,7);
% init_lens0 = fib_data0(:,9);
% fib_stress0 = fib_data0(:,12);
% fib_stretch0 = fib_data0(:,10);

stretch= load('stretch.txt');

xstress= load('xstress.txt');
ystress= load('ystress.txt');
zstress= load('zstress.txt');
xystress= load('xzstress.txt');
% 
% lens0 = init_lens0(fib_type==1);
% rads0=fib_rads0(fib_type==3);
% rads=fib_rads(fib_type==3);
% lens = init_lens(fib_type==1);
% rads01=fib_rads0(fib_type==1);
% rads1=fib_rads(fib_type==1);
% lens03 = init_lens0(fib_type==3); 
% lens3 = init_lens(fib_type==3);
% 
% origin_stress1 = fib_stress0(fib_type==1);
% stress1 = fib_stress(fib_type==1);
% origin_stress2 = fib_stress0(fib_type==3);
% stress2 = fib_stress(fib_type==3);
% 
% origin_stretch1 = fib_stretch0(fib_type==1);
% stretch1 = fib_stretch(fib_type==1);
% origin_stretch2 = fib_stretch0(fib_type==3);
% stretch2 = fib_stretch(fib_type==3);
% 
% figure
% histogram(origin_stretch1)
% hold on
% histogram(stretch1)
% xlabel('Actin Fiber Stretch')
% 
% figure
% histogram(origin_stretch2)
% hold on
% histogram(stretch2)
% xlabel('Collagen Fiber Stretch')
% 
% figure
% histogram(origin_stress1)
% hold on
% histogram(stress1)
% xlabel('Actin Fiber Stress')
% 
% figure
% histogram(origin_stress2)
% hold on
% histogram(stress2)
% xlabel('Collagen Fiber Stress')
% 
% figure
% histogram(rads0)
% hold on
% histogram(rads)
% xlabel('Collagen Fiber Radius')
% 
% figure
% histogram(rads01)
% hold on
% histogram(rads1)
% xlabel('Actin Fiber Radius')
% 
% figure
% histogram(lens0)
% hold on
% histogram(lens)
% xlabel('Actin Fiber Length')
% 
% figure
% histogram(lens03)
% hold on
% histogram(lens3)
% xlabel('Collagen Fiber Length')

figure
plot(time,xstress(:,1),'b-','LineWidth',2)
hold on
plot(time,xstress(:,2),'b--','LineWidth',2)
plot(time,xstress(:,3),'b:','LineWidth',2)
plot(time,xstress(:,4),'b-.','LineWidth',2)
plot(time,ystress(:,1),'g-','LineWidth',2)
plot(time,ystress(:,2),'g--','LineWidth',2)
plot(time,ystress(:,3),'g:','LineWidth',2)
plot(time,ystress(:,4),'g-.','LineWidth',2)
plot(time,zstress(:,1),'r-','LineWidth',2)
plot(time,zstress(:,2),'r--','LineWidth',2)
plot(time,zstress(:,3),'r:','LineWidth',2)
plot(time,zstress(:,4),'r-.','LineWidth',2)
% plot(time,xystress(:,1),'m-','LineWidth',2)
% plot(time,xystress(:,2),'m--','LineWidth',2)
% plot(time,xystress(:,3),'m:','LineWidth',2)
% plot(time,xystress(:,3),'m-.','LineWidth',2)
legend('Total','Actin','Elastin','Collagen')
xlabel('Time (h)')
ylabel('Stress [Pa]')
set(gca,'FontSize',16)

vol = stretch(:,1).*stretch(:,2).*stretch(:,3);
%vol = vol./vol(480);
figure
plot(time,stretch(:,1),'b-','LineWidth',2)
hold on
plot(time,stretch(:,2),'g--','LineWidth',2)
plot(time,stretch(:,3),'r:','LineWidth',2)
plot(time, vol, 'k-', 'lineWidth',2)
legend('1-1 Stretch','2-2 Stretch','3-3 Stretch','V/V_0')
xlabel('Time (h)')
ylabel('Stretch, \lambda')
set(gca,'FontSize',16)

figure
plot(time,phi(:,1),'b-','LineWidth',2)
hold on
plot(time,phi(:,2),'g--','LineWidth',2)
plot(time,phi(:,3),'r:','LineWidth',2)

legend('Actin','Elastin','Collagen')
xlabel('Time (h)')
ylabel('Fiber Volume Fraction, \phi')
set(gca,'FontSize',16)

figure
plot(time,orient(:,1),'b-','LineWidth',2)
hold on
plot(time,orient(:,2),'g--','LineWidth',2)
plot(time,orient(:,3),'r:','LineWidth',2)
legend('1-1','2-2','3-3')
xlabel('Time (h)')
ylabel('Orientation')
set(gca,'FontSize',16)

figure
plot(time,WSS,'k-','LineWidth',2)
xlabel('Time (h)')
ylabel('Wall Shear Stress, [Pa]')
set(gca,'FontSize',16)

% rve_stretch0 = stretch(1,:);
% rve_stretch = stretch(end,:);
% 
% plot_net_single_fib_type(nodes0, bnd_nodes0, fibers0, fib_type, rve_stretch0)
% 
% plot_net_single_fib_type(nodes, bnd_nodes, fibers, fib_type, rve_stretch)

% for p=1:101
%     E1=175e3;
%     E2=1e6;
%     A0 = pi()*(2.5e-9)^2;
%     
%     lam(p) = 0.01*(p-1) + 1;
%     
%     GS = 0.5*(lam(p)^2-1);
%     
%     sig1(p) = E1*GS;
%     sig2(p) = E2*GS;
% end
% 
% figure
% plot(lam,sig1)
% hold on
% plot(lam,sig2)
% xlabel('Stretch')
% ylabel('Stress')