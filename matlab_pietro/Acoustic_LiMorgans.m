clear
close all

load 'Meanflow_Tsin_1600_800.mat'
clear alpha 
alpha = -beta;
clear beta P_inlet x1 M1 P1 M_inlet i j

%% COMPUTING Tud
integral_p=zeros(1,N+1);
integral_m=zeros(1,N+1);
P1_p=zeros(1,N+1);
P1_m=zeros(1,N+1);
u1_p=zeros(1,N+1);
u1_m=zeros(1,N+1);

He = 1.5;

omega=He*2*pi*c(1)/L;
k0=omega/c(1);

Z1 = -1;
%u0 = u(1);
u0=1;
P0 = rho(1)*c(1)*Z1*u0;

for n=1:1:N+1
    k0=omega/c(n);
    integral_p(n)=exp(sum(-1i*omega*dx./(c(1:n)+u(1:n))));
    P1_p(n)=(rho(n)/rho(1))^0.25*(1+M(1))/(1+M(n))*(exp(gamma*M(1)-gamma/4*M(1)^2-(gamma^2-1)/3*M(1)^3))/(exp(gamma*M(n)-gamma/4*M(n)^2-(gamma^2-1)/3*M(n)^3))*integral_p(n);
    u1_p(n)=((1i*k0-(1+2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))*P1_p(n)/(rho(n)*c(n));
    
    integral_m(n)=exp(sum(1i*omega*dx./(c(1:n)-u(1:n))));
    P1_m(n)=(rho(n)/rho(1))^0.25*(1-M(1))/(1-M(n))*(exp(gamma*M(n)+gamma/4*M(n)^2-(gamma^2-1)/3*M(n)^3))/(exp(gamma*M(1)+gamma/4*M(1)^2-(gamma^2-1)/3*M(1)^3))*integral_m(n);
    u1_m(n)=(-(1i*k0+(1-2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))*P1_m(n)/(rho(n)*c(n));
end

n = 1;
k0=omega/c(1);
C1=inv([1 1; ((1i*k0-(1+2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))/(rho(n)*c(n)) (-(1i*k0+(1-2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))/(rho(n)*c(n))])*[P0; u0];

for n=1:1:N+1
    X=[P1_p(n) P1_m(n); u1_p(n) u1_m(n)]*C1;
    P1(n)=X(1);
    u1(n)=X(2);
    clear X
end

Fp = P1./(rho(1)*c(1)*u0);
Fu = u1./u0;

subplot(2,2,1)
plot(x/L,abs(Fp),'o')
ylabel('$|F_p|$','interpreter','latex')
grid on
ylim([0.95 1.2])

subplot(2,2,3)
plot(x/L,unwrap(angle(Fp))/pi,'o')
xlabel('$\frac{x}{L}$','interpreter','latex')
ylabel('$\frac{\angle{F_p}}{\pi}$','interpreter','latex')
grid on
ylim([-1.1 1.1])

subplot(2,2,2)
plot(x/L,abs(Fu),'o')
ylabel('$|F_u|$','interpreter','latex')
grid on
ylim([0.65 1.1])

subplot(2,2,4)
plot(x/L,unwrap(angle(Fu))/pi,'o')
xlabel('$\frac{x}{L}$','interpreter','latex')
ylabel('$\frac{\angle{F_u}}{\pi}$','interpreter','latex')
grid on
ylim([-1.1 1.1])

%save('Tlin1600-800.mat', '','','')
save('Tsin1600-800')