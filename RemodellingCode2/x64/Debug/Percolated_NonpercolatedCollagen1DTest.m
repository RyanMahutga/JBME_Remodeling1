%percolation test
close all
clear all

ke=1e6*pi()*10e-9^2;
kcp=7e6*pi()*100e-9^2;

kcnp=70e6*pi()*100e-9^2;

lamL=1.15;

for i=1:250
    
    lam(i) = 1+(i-1)*0.001;
    
   % nonperc collagen
   
   if lam(i) < lamL
       keq_np = ke/2;
       F_np(i) = keq_np*(lam(i)-1);
   else
       if lam(i) == lamL
           F_tr = keq_np*(lamL-1);
       end
       keq_np = 1/(1/ke+1/(ke+kcnp));
       F_np(i) = keq_np*(lam(i)-lamL)+F_tr;
   end
   
   if lam(i) < lamL
       keq_p = ke/2;
       F_p(i) = keq_p*(lam(i)-1);
   else
       keq_p = ke/2 + kcp/2;
       F_p(i) = keq_p*(lam(i)-lamL)+F_tr;
   end
   
end

figure
plot(lam,F_p.*1e9,'k--','LineWidth',2)
hold on
plot(lam,F_np.*1e9,'k:','LineWidth',2)
legend('Percolating', 'Non-Percolating')
xlabel('Stretch, \lambda')
ylabel('Force, [nN]')
set(gca,'FontSize',16)
