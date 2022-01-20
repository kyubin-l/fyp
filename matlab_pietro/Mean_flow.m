
clear
close all

%% GAS PROPERTIES
gamma=1.4;
R=287;

%% DEFINE TUBE LENGTH AND DISCRETISATION
L=1; %tube length
N=2000; %intervals along the length
x=linspace(0,L,N+1);
dx=x(2)-x(1);

%% DEFINE INITIAL CONDITIONS
M_inlet = 0.2;
P_inlet = 100000;
T1 = 1600;
T2 = 800;

%% DEFINE AREA AND TEMPERATURE DISTRIBUTION
A=zeros(1,N+1);
alpha=zeros(1,N+1);
T=zeros(1,N+1);
beta=zeros(1,N+1);

%area distribution
A=-0.*(x./L).^2+4;
%temperature distribution
T=(T2-T1).*x/L+T1;%/2*sin(5*pi*(x-x(1))/(4*L)+pi/4)+(T1+T2)/2;
%T = (T1-T2)/2*sin(5*pi*(x-x(1))/(4*L)+pi/4)+(T1+T2)/2;

%% COMPUTING ALPHA AND BETA DISTRIBUTIONS
for i=2:length(A)-1
    alpha(i)=(A(i+1)-A(i-1))/(2*A(i)*dx);
    beta(i)=(T(i+1)-T(i-1))/(2*T(i)*dx);   
end
alpha(1)=(A(2)-A(1))/(A(1)*dx);
beta(1)=(T(2)-T(1))/(T(1)*dx);
alpha(end)=(A(end)-A(end-1))/(A(end)*dx);
beta(end)=(T(end)-T(end-1))/(T(end)*dx);

%% COMPUTING MACH DISTRIBUTION
M(1)=M_inlet;
for j=1:N
    [x1,M1] = ode45(@(x1,M1) M1*(beta(j)*(1+gamma*M1^2)-2*alpha(j))/(2*(1-gamma*M1^2)),[x(j) x(j+1)],M(j));
    M(j+1)=M1(end);
end

%% COMPUTING MEAN FLOW DISTRIBUTIONS
c=sqrt(gamma*R*T); %speed of sound
u=M.*c;
P(1)=P_inlet;
for j=1:N
    [x1,P1] = ode45(@(x1,P1) -P1*(gamma*M(j)^2*(beta(j)-alpha(j)))/(1-gamma*M(j)^2),[x(j) x(j+1)],P(j));
    P(j+1)=P1(end);
end
rho=P./(R.*T);

%% PLOT MEAN FLOW PROPERTIES
figure
sgtitle('Mean flow properties')
subplot(4,1,1)
plot(x,M)
ylabel('Mach')
grid on
subplot(4,1,2)
plot(x,u)
ylabel('Velocity [m/s]')
grid on
subplot(4,1,3)
plot(x,P)
ylabel('Pressure [Pa]')
grid on
subplot(4,1,4)
plot(x,rho)
ylabel('Density [kg/m3]')
grid on
xlabel('Distance from inlet [m]')

save('Meanflow_Tlin_1600_800')
