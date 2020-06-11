close all

R=8.314;
T=310;

c_star=150;

k = 200;
k2=381;
k3=0.241;

for p=1:101

phim(p) =  (p-1)/100;   
    
c_fcd(p) = k * phim(p);

pressure(p) = R * T * (sqrt(c_fcd(p)*c_fcd(p) + 4 * c_star*c_star) - 2 * c_star) + k2*phim(p)*c_fcd(p) + k3*phim(p)^2*c_fcd(p)^2; 

end

figure
plot(phim,pressure./1000,'k-','LineWidth',2)
xlabel('Phi')
ylabel('Pressure [kPa]')