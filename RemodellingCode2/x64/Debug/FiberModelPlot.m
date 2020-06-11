% Fiber model
% close all
% clear all
% clc

fib_mod = 0.7e9;
fib_area = pi()*100e-9*100e-9;

R0 = 5.8; % [nm] roughly radius of collagen microfibril
r0 = 1.6; % [nm] roughly a collagen triple helix
H0 = 67.4; %[nm] d-pattern banding in collagen molecules
PI=pi();

L0 = sqrt((2 * PI*R0)*(2 * PI*R0) + H0*H0);

lambda_bar = L0 / H0;

E_bar = fib_mod * H0*H0 / (L0*(H0 + (1 + 37 / (6 * PI*PI) + 2 * (L0*L0) / ((PI*r0)*(PI*r0)))*(L0 - H0)));

sigma = 0;

for k=1:1501
    lambda = 0.65+(0.001*(k-1));
    
    if (lambda <= lambda_bar) % see above paper for more information
        H = lambda * H0;
        R = sqrt(L0*L0 - H*H) / (2 * PI);
        eta = (R*R + H*H) / (L0*H*(1 + 4 * R*R / (r0*r0) + 6 * (20 / 9 + R*R / (r0*r0))*R*R / (H*H)));
        
        dHdlam = H0;
        
        dRdH = -H/(2*PI*sqrt((L0*L0 - H*H)));
        
        dEtadH = -(H*H + R*R) / ((H*H) * L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) +...
            2 / (L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) +...
            12 * (R*R) * ((H*H) + (R*R))*((R*R) / (r0*r0) + 20 / 9) / ((H*H*H*H) * L0*((6 *(R*R) *...
            ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)*(6 * (R*R) * ((R*R) / ...
            (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)));
        
        dEtadR = 2 * R / (H*L0*(6 * (R*R) * ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)) -...
            ((H*H) + (R*R))*(12 * (R*R) / ((r0*r0) * (H*H)) + 12 * R*((R*R) /...
            (r0*r0) + 20 / 9) / (H*H) + 8 * R / (r0*r0)) / (H*L0*((6 * (R*R) *...
            ((R*R) / (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)*(6 * (R*R) * ((R*R) /...
            (r0*r0) + 20 / 9) / (H*H) + 4 * (R*R) / (r0*r0) + 1)));
        
        dEtadlam = dEtadH * dHdlam + dEtadR * dRdH*dHdlam;
        
        sigma = eta * E_bar*(lambda - 1);
        dsigma = (eta*E_bar + E_bar * (lambda - 1)*dEtadlam);
        dFdlam = fib_area * (eta*E_bar + E_bar * (lambda - 1)*dEtadlam);
        
    else
        
        sigma = E_bar * (lambda_bar - 1) + fib_mod * (lambda / lambda_bar - 1);
        dsigma = fib_mod / lambda_bar;
        dFdlam = fib_mod * fib_area / lambda_bar;
    end
    F(k) = sigma*fib_area;
    stress(k) = sigma;
    dstress(k) = dsigma;
    dF(k) = dFdlam;
    lam(k)=lambda;
    
    GS = (lambda-1);
    
    Ee =  1e6;
    if lambda < 1
        Ee = 1e-2*Ee ;
    end
    
    stressE(k) = Ee*GS;
    Fe(k) = stressE(k)*fib_area;
    
    Em = 3000e3;
    Em2=2000e3;
    if lambda < 1
        Em = 1e-2*Em ;
        Em2=1e-2*Em2;
    end
    lam_max = 1;
    lam_min = 0.85;
    S=500e3;
    S2=200e3;
    
    stressM_p(k) = Em*GS;
    Fmp(k) = Em*GS*fib_area/100;
    f_lam = 1-(lam_max-lambda)^2/(lam_max-lam_min)^2;
        
    

    stressM_a(k) = S*f_lam;
    Fma(k) = S*f_lam*fib_area/100;
    if stressM_a(k)<0
           stressM_a(k) = 0;
    Fma(k) = 0; 
    end
    stressM(k) = stressM_a(k)+stressM_p(k);
    Fm(k) = Fma(k);
    
    stressM2_p(k) = Em2*GS;
    Fmp2(k) = Em2*GS*fib_area/100;
    f_lam2 = 1-(lam_max-lambda)^2/(lam_max-lam_min)^2;
%         
    stressM2_a(k) = S2*f_lam2;
    Fma2(k) = S2*f_lam2*fib_area/100;
    stressM2(k) = stressM2_a(k)+stressM2_p(k);
    Fm2(k) = Fma2(k);
end

idx =find(stress > 50e3,1,'first');
idx2 = find(stressM > 450e3,1,'first');
% 
% figure
% plot(lam,F.*1e9,'b-','LineWidth', 2)
% hold on
% plot(lam,Fe.*1e9,'g-','LineWidth', 2)
% plot(lam,Fm.*1e9,'r-','LineWidth', 2)
% plot(lam,Fma.*1e9,'r:','LineWidth', 2)
% % plot(lam,Fmp.*1e9,'r--','LineWidth', 2)
% plot(lam,Fm2.*1e9,'m-','LineWidth', 2)
% plot(lam,Fma2.*1e9,'m:','LineWidth', 2)
% % plot(lam,Fmp2.*1e9,'m--','LineWidth', 2)
% xlabel('Stretch')
% ylabel('Fiber Force (nN)')
% set(gca,'FontSize',18)
% legend('Collagen','Elastin','Actin','Active','Passive')
% axis([0.65 1.55 -0.5 10.5])

% figure
% plot(lam, dF.*1e9,'r--','LineWidth', 2)
% xlabel('Stretch')
% ylabel('Tangent Modulus (nN)')
% set(gca,'FontSize',18)
% 
% comp_modulus = mean(dF(1:50));
% tens_modulus = mean(dF(52:101));
% 
% modulus_ratio = tens_modulus/comp_modulus;

figure
plot(lam, stress./1e3,'-','Color',[0 0 0.8],'LineWidth', 2)
hold on
plot(lam,stressE./1e3,'-','Color',[0 0.7 0],'LineWidth',2)
plot(lam,stressM./1e3,'-','Color',[0.8 0 0],'LineWidth',2)
plot(lam,stressM_a./1e3,':','Color',[0.8 0 0],'LineWidth',1)
plot(lam,stressM_p./1e3,'--','Color',[0.8 0 0],'LineWidth',1)
% plot(lam,stressM2./1e3,'m-','LineWidth',2)
% plot(lam,stressM2_a./1e3,'m:','LineWidth',1)
% plot(lam,stressM2_p./1e3,'m--','LineWidth',1)
% plot(lam(idx), stress(idx)/1e9,'x','Color',[0 0 0.8],'MarkerSize',10,'LineWidth',2)
% plot(lam(idx2),stressM(idx2)./1e9,'x','Color',[0.8 0 0],'MarkerSize',10,'LineWidth',2)
%plot([1.1 1.1],[-10 10],'k:','LineWidth',2)
%plot([1.13 1.13],[-10 10],'k:','LineWidth',2)
%  plot(1.167, 0.2,'ro','LineWidth',2,'MarkerSize', 8);
%  plot(1.1375, 1.75,'ro','LineWidth',2,'MarkerSize', 8);
%  plot(1.167, 0.2,'rx','LineWidth',2,'MarkerSize', 6);
%  plot(1.1375, 1.75,'rx','LineWidth',2,'MarkerSize', 6);
xlabel('Stretch')
ylabel('Fiber Stress [MPa]')
legend('Collagen','Elastin','Actin','Active','Passive')
set(gca,'FontSize',18)
axis([0.8 1.3 -50 1500])

%plot(, 2e6,'ro','MarkerSize', 3);

% figure
% plot(lam, dstress./1e6,'r--','LineWidth', 2)
% xlabel('Stretch')
% ylabel('Tangent Modulus (MPa)')
% set(gca,'FontSize',18)


