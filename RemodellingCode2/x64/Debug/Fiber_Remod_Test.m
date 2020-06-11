clear all
close all
clc

rad = [2e-9,20e-9, 2e-9];

init_len = [1,1,1];

G1=1e-8;
G2=5e6;
G3=0.85e-8;
G4=3.5e6;

tau=[1440,60e3];
sig_targ=[50e3,35e3];

E=[500e6,5e6, 300e6];

stress(1,:)=[0,0,0];

fib_rads = rad;
init_lens = init_len;

Ftarget = 1e-9;
A(1)=1;

k_r=1;

time(1)=0;
err(1)=0;
for t = 1:4000
    area = pi().*rad.^2;
    
    L = fminsearch(@(L) force_equilibrium(L,init_len,area,Ftarget,E), 1.1);
    err(t+1) = force_equilibrium(L,init_len,area,Ftarget,E);
    
    if t==1
        V = L*A(1);
        phi = init_len(3)/V;
    end
    
    A(t+1) = init_len(3)/(L*phi);
    
    if t<100
        dt=0.001;
    elseif t<200
        dt=0.1;
    elseif t<300
        dt=0.5;
    elseif t<2000
        dt=1;
    elseif t<2100
        dt=0.001;
        Ftarget=1.5e-9;
    elseif t<2200
        dt=0.1;
    elseif t<2300
        dt=0.5;
    else
        dt=1;
    end
  
    time(t+1) = time(t) + dt;
    
    lambda = L./init_len;
    sigma = E.*(lambda-1);
    
    for k=1:length(init_len)     
        if k==1 % collagen
            if sigma(k) >0
                Kdep = G3*exp(sigma(k)/G4);
                Kdeg = G1*exp(-sigma(k)/G2);
                
                dL = rad(k)*init_len(k)/(2*k_r*init_len(k))*(Kdep-Kdeg); %1/tau(1)*(sigma(k)/sig_targ(1)-1)*rad(k)*dt;
                dR= k_r*dL;
            else
                dL = rad(k)*init_len(k)/(2*k_r*init_len(k))*(G3-G1);
                dR = k_r*dL ; %-1/tau(1)*rad(k)*dt;
            end
            init_len(k) = init_len(k) + dL ; %*(rad(k)+dR)^2/(rad(k)^2+(2*rad(k)*dR+dR^2)/lambda(k));            
            rad(k)=rad(k)+dR;
        elseif k==3 % actin
            if sigma(k) > 0
                dL = 0;%1/tau(2)*(sigma(k)/sig_targ(2)-1)*init_len(k)*dt;
            else
                dL = 0;%-1/tau(2)*init_len(k)*dt;
            end
            init_len(k) = init_len(k) + dL ; %*(rad(k)+dR)^2/(rad(k)^2+(2*rad(k)*dR+dR^2)/lambda(k));            
            %rad(k)=rad(k)+dR;
        end
    end
    stress(t+1,:)=sigma;
    fib_rads(t+1,:) = rad;
    init_lens(t+1,:) = init_len;
    
end

figure
plot(time,fib_rads)
legend('Collagen','Elastin','Actin')
xlabel('Time')
ylabel('Radii')

figure
plot(time,init_lens)
legend('Collagen','Elastin','Actin')
xlabel('Time')
ylabel('Lens')

figure
plot(time,stress)
legend('Collagen','Elastin','Actin')
xlabel('Time')
ylabel('Stress')

figure
plot(time,err)


function [err] = force_equilibrium(L, init_len, area, Ftarget,E)

lambda = L./init_len;

F=0;
for k=1:length(init_len)
    if k==1
        if lambda > 1.14
            sigma = E(k)*(lambda(k)-1);
        else
            sigma = 1e-12*E(k)*(lambda(k)-1);
        end
    else
        if lambda > 1
            sigma = E(k)*(lambda(k)-1);
        else
            sigma = 1e-12*E(k)*(lambda(k)-1);
        end
    end
    F = F + sigma*area(k);
end

err = (F-Ftarget)^2;
end