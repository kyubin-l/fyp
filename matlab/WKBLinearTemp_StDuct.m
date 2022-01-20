%% ACOUSTIC FIELD FOR AN NON-ISENTROPIC FLOW THROUGH A QUASI 1D DUC 
clc
clear all
close all
format long 
tic
%sigma_hat_in, M_in, ro, ri, L, He, T_in, T_end,
scale = 0.0001;
% Scales of solvers (all same at the moment)
mean_scale = scale; Beta_scale = scale; % 1 step per mm

% Geometric description of the duct 
L = 1;                  % Length of the duct in m
xnum = (0:mean_scale:L)';
 
% Temperature profile
T_in = 1600; T_end = 800; % Temperature drop of 800K/m. Here the length is 0.5m

% Other inlet conditions
Z = -1;           % Impedance condition
M_in = 0.2;      % Inlet Mach number, varies with space
He = 1.5;        % Helmholtz number based on length

% Flow constants
Gamma = 1.4;             % Ratio of specific heats for air, treated as a constant
R_g = 287;                % Gas constant for air in Joules/kg.Kelvin

% Mean flow at the input            
P_in = 1 * 10^5;         % Pressure at inlet in Pascals
Rho_in = P_in/(R_g*T_in);        
U_in = M_in * sqrt(Gamma * R_g * T_in); % Inlet Velocity of the flow in m/s

% Pertubation inputs
p_hat_in = 100;                  % Inlet acoustic disturbance (in Pa)                         
Omega = He*2*pi*sqrt(Gamma * R_g * T_in)/L;     % 2pi * Frequency
Freq = Omega/(2*pi);

%% SOLVING THE STEADY CONSERVATION EQUATIONS - Mean flow
% y(2) - P; y(1) - U
xspan = (0:mean_scale:L)'; y0 = [U_in P_in Rho_in];
opts = odeset('RelTol',1e-15,'AbsTol',1e-20);

% Temperature Distribution
T_x = (0:scale:L)';
T(1,1) = T_in;
Beta_span = (0:Beta_scale:L)';
% Linear
% Beta(1,1) = (T_end - T(1))/(L*T(1));
% Sinusoidal
Beta(1,1) = 1/T(1)*((T(1) - T_end)/2)*cos(5*pi*xspan(1)/(4*L) + pi/4) * 5*pi/(4*L);

for i = 1:size(xspan,1)
    % Linear dist considered
    % T(i,1) = T_in + (T_end - T_in)*(xspan(i)/L);
    % Beta(i,1) =  (1/T(i,1))*(T_end - T_in)/L;
    
    % Sinusoidal
    T(i,1) = ((T_in - T_end)/2)*sin(5*pi*xspan(i)/(4*L) + pi/4) + (T_in + T_end)/2;
    Beta(i,1) = 1/T(i,1)*(T_in - T_end)/2*cos(5*pi*xspan(i)/(4*L) + pi/4) * 5*pi/(4*L);
end
% Solving
[x,y] = ode45(@(x,y) odexMeanIsothermalfcn(x,y,R_g,Gamma,T,T_x,Beta,Beta_span), xspan, y0, opts);
% Exporting Results
U = y(:,1);
P = y(:,2);
Rho = y(:,3);
y(:,4) = T;
M = U./sqrt(Gamma*R_g.*T);
c = sqrt(Gamma .* R_g .* T);
ko = Omega./c;           % Wave number

%% SOLVING 2 LINEARISED EULER EQUATIONS
% U^ - y(1); P^ - y(2)   
u_hat_in = p_hat_in/(Rho(1)*c(1)*Z);       % Velocity perturbation in Pa
xnewspan = (0:scale:L)';                         % Where to evaluate the final solution   
ynew0 = [u_hat_in p_hat_in];
optsnew = odeset('RelTol',1e-15,'AbsTol',1e-20);

% For interpolation in ODE45 routine
M_x = (0:scale:L)';
U_x = (0:scale:L)';   
P_x = (0:scale:L)';
T_x = (0:scale:L)'; 
Rho_x = (0:scale:L)';
Beta_x = Beta_span;

[xnew,ynew] = ode45(@(xnew,ynew) odex2LEEfcnNonIsentropic(xnew,ynew,M_x,U_x,P_x,T_x,Rho_x,M,U,P,T,R_g,Rho,Gamma,Omega,Beta,Beta_x), xnewspan, ynew0, optsnew);
U_hat_2LEE = ynew(:,1);                                        
P_hat_2LEE = ynew(:,2);
Fp_2LEE = P_hat_2LEE./(Rho(1) * c(1) * u_hat_in); % Normalising
Fu_2LEE = U_hat_2LEE./u_hat_in;                       % Normalising

Fp_mag_2LEE = abs(Fp_2LEE);       % Magnitude of Fp  
Fu_mag_2LEE = abs(Fu_2LEE);       % Magnitude of Fu

Fp_phase_2LEE = angle(Fp_2LEE)./pi;   % Phase of Fp    
Fu_phase_2LEE = angle(Fu_2LEE)./pi;   % Phase of Fu  
%% Semi Analytical Non Isentropic definition
Alpha = (1./Rho) .* gradient(Rho)./gradient(x);
for i = 1:size(xnew,1)
       
    if i==1
    A_plus(i,1) =  (1/(1i*ko(i,1) - (Alpha(i,1))*M(i,1))) * (1i*ko(i,1) - (Alpha(i,1)/4)*(1 + 2*(1+Gamma)*M(i,1)+ (3*Gamma-7)*M(i,1)^2));     
    A_minus(i,1) =  (1/(1i*ko(i,1) - (Alpha(i,1))*M(i,1))) * (1i*ko(i,1) + (Alpha(i,1)/4)*(1 - 2*(1+Gamma)*M(i,1)+ (3*Gamma-7)*M(i,1)^2));
    
    C1_plus_Li = p_hat_in * (A_minus(1,1) + (1/Z))/(A_plus(1,1) + A_minus(1,1));
    C1_minus_Li = p_hat_in * (A_plus(1,1) - (1/Z))/(A_plus(1,1) + A_minus(1,1));
    
    P_hat_Analytical_Li(1,1) = p_hat_in;
    U_hat_Analytical_Li(1,1) = u_hat_in;
    P1_plus_Li(1,1) = 1; P1_minus_Li(1,1) = 1;
    else
 
    A_plus(i,1) =  (1/(1i*ko(i,1) - (Alpha(i,1))*M(i,1))) * (1i*ko(i,1) - (Alpha(i,1)/4)*(1 + 2*(1+Gamma)*M(i,1)+ (3*Gamma-7)*M(i,1)^2));     
    A_minus(i,1) =  (1/(1i*ko(i,1) - (Alpha(i,1))*M(i,1))) * (1i*ko(i,1) + (Alpha(i,1)/4)*(1 - 2*(1+Gamma)*M(i,1)+ (3*Gamma-7)*M(i,1)^2));
    
    I1 = 1./(c(1:i,1) + U(1:i,1)); I2 = 1./(c(1:i,1) - U(1:i,1));
    
    Phase_P1_plus_Li = exp(-1i*Omega * scale * cumsum(I1(1:i,1)));
    Phase_P1_minus_Li = exp(1i*Omega * scale * cumsum(I2(1:i,1)));
    
    Phase_P1_plus_Li = exp(-1i*Omega * scale * cumsum(I1(1:i,1)));
    Phase_P1_minus_Li = exp(1i*Omega * scale * cumsum(I2(1:i,1)));
    
    Mag_P1_plus(i,1) = ((Rho(i,1)/Rho(1,1))^0.25) * ((1+M(1,1))/(1+M(i,1))) * (exp( Gamma*(M(1,1)-M(i,1)) - (Gamma/4)*(M(1,1)^2 - M(i,1)^2) - ((Gamma^2 -1)/3) * (M(1,1)^3 - M(i,1)^3  )));    
    Mag_P1_minus(i,1) = ((Rho(i,1)/Rho(1,1))^0.25) * ((1-M(1,1))/(1-M(i,1))) *(exp( Gamma*(M(i,1)-M(1,1)) + (Gamma/4)*(M(i,1)^2 - M(1,1)^2) - ((Gamma^2 -1)/3) * (M(i,1)^3 - M(1,1)^3  )));        
    
    P1_plus_Li(i,1) = Mag_P1_plus(i,1) * Phase_P1_plus_Li(i,1);
    P1_minus_Li(i,1) = Mag_P1_minus(i,1) * Phase_P1_minus_Li(i,1);

    P_hat_Analytical_Li(i,1) = C1_plus_Li * P1_plus_Li(i,1) + C1_minus_Li * P1_minus_Li(i,1);
    U_hat_Analytical_Li(i,1) = (1/(Rho(i,1) * c(i,1))) * (A_plus(i,1) * C1_plus_Li* P1_plus_Li(i,1) - A_minus(i,1) * C1_minus_Li * P1_minus_Li(i,1));
    end  
end
Fp_Analytical = P_hat_Analytical_Li./(Rho(1,1) * c(1,1) * u_hat_in);
Fu_Analytical = U_hat_Analytical_Li./(u_hat_in);

Fp_mag_Analytical = abs(Fp_Analytical);           % Magnitude of Fp  
Fu_mag_Analytical = abs(Fu_Analytical);           % Magnitude of Fu
Fp_phase_Analytical = angle(Fp_Analytical)./pi;   % Phase of Fp    
Fu_phase_Analytical = angle(Fu_Analytical)./pi;   % Phase of Fu 

%% Saving files and data

save('WKB_M0.2_Tsin.mat', 'x', 'Fp_2LEE', 'Fu_2LEE', 'Fp_Analytical', 'Fu_Analytical')

%% Error Definition
% For Analytical technique
EpsilonP_ANA = sqrt(sum((abs(Fp_2LEE - Fp_Analytical)).^2)./sum((abs(Fp_2LEE)).^2));
EpsilonU_ANA = sqrt(sum((abs(Fu_2LEE - Fu_Analytical)).^2)./sum((abs(Fu_2LEE)).^2));
%% Plotting the results
ylabels{1}='Velocity in m/s';
ylabels{2}='Pressure in Pa';
ylabels{3}='Temperature in K';
xlabel('x Location');
[ax,hlines] = multiplotyyy({x,y(:,3)},{x,y(:,2)},{x,y(:,1)},ylabels);
legend(cat(1,hlines{:}),'Velocity in m/s','Pressure','Temperature in K','location','e')
title('Mean flow solution');
ax(2).Position = ax(1).Position;
ax(3).Position = ax(1).Position;
offset = 0.7/5.5;
ax(3).Position(3) = ax(3).Position(3) + offset;
ax(3).XLim = [ ax(1).XLim(1) ax(1).XLim(2)*(1 + 1/(ax(1).Position(3)) * offset)];
grid on; 
hold off

% Plotting normalised perturbation field of pressure and velocity, as phase and gain
b = ceil(length(xnew)/20);
figure;
subplot(2,2,1)
hold on
plot(xnew/L,Fp_mag_2LEE,'r-*','MarkerIndices',round(b/3):b:length(xnew))
plot(xnew/L,Fp_mag_Analytical,'b:^','MarkerIndices',2*round(b/3):b:length(xnew))
ylim([0.95 1.2])
legend('2LEE','Model')
xlabel('x Location');
ylabel('F_p Gain'); % ESTABLISHING THE MEAN PROPERTIES
grid on; 
hold off
subplot(2,2,2)
hold on
plot(xnew/L,Fp_phase_2LEE,'r-*','MarkerIndices',round(b/3):b:length(xnew))
plot(xnew/L,Fp_phase_Analytical,'b:^','MarkerIndices',2*round(b/3):b:length(xnew))

xlim([0,1])
xlabel('x Location');
ylabel('F_p Phase'); 
grid on; 
hold off

subplot(2,2,3)
hold on
plot(xnew/L,Fu_mag_2LEE,'r-*','MarkerIndices',round(b/3):b:length(xnew))
plot(xnew/L,Fu_mag_Analytical,'b:^','MarkerIndices',2*round(b/3):b:length(xnew))
xlim([0,1])
xlabel('x Location');suptitle(['Acoustic solution for Mach number of ',num2str(M_in),' at inlet and ',num2str(ceil(100*M(end,1))/100),'at exit'])
ylim([0.65 1.05])
ylabel('F_u Gain');

grid on; 
hold off

subplot(2,2,4)
hold on
plot(xnew/L,Fu_phase_2LEE,'r-*','MarkerIndices',round(b/3):b:length(xnew))
plot(xnew/L,Fu_phase_Analytical,'b:^','MarkerIndices',2*round(b/3):b:length(xnew))
xlabel('x Location');
xlim([0,1])
ylabel('F_u Phase');  
grid on; 

hold off
% Functions Definition
%%
function dydx = odexMeanIsothermalfcn(x,y,R_g,Gamma,T,T_x,Beta,Beta_x)
% y(2) - P; y(1) - U
T = interp1(T_x,T,x);
Beta = interp1(Beta_x,Beta,x);
dydx = zeros(3,1);
dydx(1) = (Beta)*y(1)*Gamma*R_g*T/(Gamma*R_g*T - (Gamma* y(1)^2));           % dU/dx = (Beta) U /(Gam R T - Gam U^2)
dydx(2) = ((- Beta)* Gamma* y(2) * y(1)^2)/(Gamma*R_g*T - Gamma*(y(1)^2));   % dpdx = alpha Gam P U^2/(Gam R T1 - Gam U^2)
dydx(3) = -(y(2)/(R_g*T)) * (Beta*Gamma*R_g*T)/(Gamma*R_g*T - Gamma*y(1)^2);
end

function dyprimedx = odex2LEEfcnNonIsentropic(xnew,ynew,M_x,U_x,P_x,T_x,Rho_x,M,U,P,T,R_g,Rho,Gamma,Omega,Beta,Beta_x)  
M = interp1(M_x,M,xnew);                           % Interpolate the data set 
P = interp1(P_x,P,xnew);
U = interp1(U_x,U,xnew);
T = interp1(T_x,T,xnew);
Rho = P/(R_g* T);
Beta = interp1(Beta_x,Beta,xnew);

% U^ - y(1); P^ - y(2) sig^ - y(3)

c = sqrt(Gamma * R_g*T); 
ko = Omega/c;
% Ap^ + Bdp^ + Cu^ + C2dU^
A = (1i*Omega + (Gamma^2)* (M^2)*U*((Beta)/(1 - Gamma*M^2)) + Beta * Gamma * U);
B = (U); 
C = (-Rho*(c^2))*((Beta)/(1 - Gamma*M^2)) * (1 + M^2 - Gamma* (M^2)) + Rho*(c^2)*Beta;
C2 = Rho*(c^2);

% Dp^ + Edp^ + Fu^ + F2dU^
D = (M^2/Rho)*((Beta)/(1-Gamma*(M^2)));       % Coefficient of p_hat
E = 1/(Rho);                                                       % Coefficient of dp_hat
F = 1i*Omega + U*((Beta)/(1-Gamma*(M^2)));    % Coefficient of u_hat
F2 = U;                                                             % Coefficient of du_hat

An = A/C2; Bn = B/C2; Cn = C/C2;
Dn = D/F2; En = E/F2; Fn = F/F2; 

Gn = (An - Dn)/(En - Bn); Hn = (Cn - Fn)/(En - Bn);   % dp_hat = Gn*p_hat + Hn*u_hat

In = A/B; Jn = C/B; Kn = C2/B;
Ln = D/E; Nn = F/E; On = F2/E; 

Qn = (Ln - In)/(Kn - On); Sn = (Nn - Jn)/(Kn - On); 
dyprimedx = zeros(2,1);
  dyprimedx(1) = Qn* ynew(2) + Sn* ynew(1);
  dyprimedx(2) = Gn* ynew(2) + Hn* ynew(1);

end















