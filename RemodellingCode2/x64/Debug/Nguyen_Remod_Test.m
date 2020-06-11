clear all
close all
clc

L_def = 1.15; % change to zero for force control
Ftarget = 2e-8;

rad = 80e-9;

init_len = 1.1;

G1=3*1.1111e-3;
G2=10e6;
G3=3*0.8333e-3;
G4=10e6;

E=700e6;

fib_rads = rad;
init_lens = init_len;

M=5;
k=0.5e-3;
c = 0.5e9;
time(1)=0;
for t = 1:6000
    area = pi()*rad^2;
    
    if L_def == 0
        L = fminsearch(@(L) force_equilibrium(L,init_len,area,Ftarget,E), 1.1);
    else
        L=L_def;
    end
    
    lambda = L/init_len;
    sigma = E*(lambda-1);
    
    sig(t)=sigma;
    
    if t<1000
        dt=0.001;
    elseif t<2000
        dt=0.01;
    elseif t<3000
        dt=0.05;
    elseif t<5000
        dt=0.1;
    elseif t<7500
        dt=0.5;
    else 
        dt=1.0;
    end
    
    time(t+1) = time(t)+dt;
    
    if sigma >0
        Kdep = G3*exp(sigma/G4);
        Kdeg = G1*exp(-sigma/G2);
        
        dR = (Kdep-Kdeg)*rad*dt;
        
    else
        dR = (G3-G1)*rad*dt;
    end
    
    if dR>0
    %dL = k*init_len*10e-6/rad*(exp(sigma/c)-1)*dt;
    init_len = (init_len * (((rad*rad/M + 2.0*rad*dR + dR * dR)) / (rad*rad/M + (2.0 * rad*dR + dR * dR) / lambda)));
    %init_len = init_len*(rad+dR)^2/(rad^2+(2*rad*dR+dR^2)/lambda);% + dL;
    end
    
    rad=rad+dR;
    
    fib_rads(t+1) = rad;
    init_lens(t+1) = init_len;
    
end

figure
plot(time,fib_rads,'k-','LineWidth',2)
xlabel('Time')
ylabel('Fiber Radius')

figure
plot(time,init_lens,'k--','LineWidth',2)
xlabel('Time')
ylabel('Fiber Length')

plot(time(2:end),sig,'k--','LineWidth', 2)
xlabel('Time')
ylabel('Stress')

function [err] = force_equilibrium(L, init_len, area, Ftarget,E)

lambda = L/init_len;

if lambda >1
    sigma = E*(lambda-1);
else
    sigma = 1e-12*E*(lambda-1);
end

F = sigma*area;

err = (F-Ftarget)^2*1e24;
end