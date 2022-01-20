clear
close all

load 'Meanflow_Tlin_1600_800.mat'
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
k0_=zeros(1,N+1);

TubeDiameter = 0.001;

He = 1;

omega=He*2*pi*c(1)/L;
k0=omega/c(1);

R = 287;
Pr = 0.71;
gamma = 1.4;

mu_ref = 1.716e-5; %kg/(ms)
T_ref = 273.15; %K
Cs = 110.4; %K

Z1 = -1;
u0 = 1;
P0 = rho(1)*c(1)*Z1*u0;

v = zeros(1,100);
v(T<=1000) = (mu_ref.*(T(T<=1000)./T_ref).^1.5.*(T_ref+Cs)./(T(T<=1000)+Cs)) ./ rho(T<=1000);
v(T>1000) = (2.653e-8.*T(T>1000)+1.573e-5) ./ rho(T>1000);

for n=1:1:N+1
    k=omega/c(n);
     
    Sh=TubeDiameter*(omega/v(n))^0.5;
    k0 = omega./(c(n)*(1+M(n))) * (1 + (1-1i)/(Sh*2^0.5) * (1 + (gamma-1)/(Pr^0.5)) - 1i/Sh^2 * (1 + (gamma-1)/(Pr^0.5) - 0.5*gamma*(gamma-1)/Pr));
    k0_(n) = omega./(c(n)*(1+M(n))) * (1 + (1-1i)/(Sh*2^0.5) * (1 + (gamma-1)/(Pr^0.5)) - 1i/Sh^2 * (1 + (gamma-1)/(Pr^0.5) - 0.5*gamma*(gamma-1)/Pr));

    integral_p(n)=exp(sum(-1i*omega*dx./(c(1:n)+u(1:n))));
    P1_p(n)=(rho(n)/rho(1))^0.25*(1+M(1))/(1+M(n))*(exp(gamma*M(1)-gamma/4*M(1)^2-(gamma^2-1)/3*M(1)^3))/(exp(gamma*M(n)-gamma/4*M(n)^2-(gamma^2-1)/3*M(n)^3))*integral_p(n);
    u1_p(n)=((1i*k0-(1+2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))*P1_p(n)/(rho(n)*c(n));
    
    integral_m(n)=exp(sum(1i*omega*dx./(c(1:n)-u(1:n))));
    P1_m(n)=(rho(n)/rho(1))^0.25*(1-M(1))/(1-M(n))*(exp(gamma*M(n)+gamma/4*M(n)^2-(gamma^2-1)/3*M(n)^3))/(exp(gamma*M(1)+gamma/4*M(1)^2-(gamma^2-1)/3*M(1)^3))*integral_m(n);
    u1_m(n)=(-(1i*k0+(1-2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))*P1_m(n)/(rho(n)*c(n));
end

n = 1;
k0 = k0_(1);
C1=inv([1 1; ((1i*k0-(1+2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))/(rho(n)*c(n)) (-(1i*k0+(1-2*(1+gamma)*M(n)+(3*gamma-7)*M(n)^2)*alpha(n)/4)/(1i*k0-alpha(n)*M(n)))/(rho(n)*c(n))])*[P0; u0];

for n=1:1:N+1
    X=[P1_p(n) P1_m(n); u1_p(n) u1_m(n)]*C1;
    P1(n)=X(1);
    u1(n)=X(2);
    clear X
end

Fp_c = P1./(rho(1)*c(1)*u0);
Fu_c = u1./u0;

%%
%%%%%%

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
plot(x,abs(Fp),x,abs(Fp_c),'o')
ylabel('$|F_p|$','interpreter','latex')
grid on
ylim([0.95 1.2])

subplot(2,2,3)
plot(x,angle(Fp)/pi,x,angle(Fp_c)/pi,'o')
xlabel('$\frac{x}{L}$','interpreter','latex')
ylabel('$\frac{\angle{F_p}}{\pi}$','interpreter','latex')
grid on
ylim([-1.1 1.1])

subplot(2,2,2)
plot(x,abs(Fu),x,abs(Fu_c),'o')
ylabel('$|F_u|$','interpreter','latex')
grid on
ylim([0.65 1.1])

subplot(2,2,4)
plot(x,angle(Fu)/pi,x,angle(Fu_c)/pi,'o')
xlabel('$\frac{x}{L}$','interpreter','latex')
ylabel('$\frac{\angle{F_u}}{\pi}$','interpreter','latex')
grid on
ylim([-1.1 1.1])

save('Tlin1600-800_damping')