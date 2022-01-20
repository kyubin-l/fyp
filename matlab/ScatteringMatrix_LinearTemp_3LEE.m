% Mean flow solver
% *Function*: MeanFlowScatteringMatrixNon-IsentropicConical
% To evaluate the mean flow and there-by the flow Mach number at the exit of the duct.
% *Inputs* : 
% b_star= (T_end - T_in) / T_in = T_star_end - 1; Since T_star = T/T_in
% Mu, - flow Machnumber at the inlet
% Output:
% Machnumber of the flow at the outlet
% Theta = Ratio of areas at inlet to outlet
% Xi = Ratio of sound speed at inlet to outlet
% Alpha_new = 1/Pressure downstream
% Beta_star = (1/T_star)*dT_star/dx* - Normalised Temp gradient
% Note ko_star_in = 0.5; % Wavenumber ko_star = Helmholtz number = omega * L/c_in

function [S_3LEE] = ScatteringMatrix_LinearTemp_3LEE(M_in, b_star,ko_star_in)
Gamma = 1.4; scale = 0.001;
xspan = (0 : scale : 1)';
[M, Beta_star] = MeanFlowForScatteringMatrixNonIsentropicLinear(b_star, M_in, scale);

T_star = (1 + b_star * xspan);
c_star = sqrt(T_star);  % c = sqrt(Gamma * R_g * T) => c* * c_in = sqrt( Gamma * R_g * T* * T_in)
% => c* = sqrt(T*)
ko = ko_star_in / c_star; % also referred as He = omega L/c
%% Full Scattering Matrix Evaluation - Non-Isentropic flow
% For interpolation in ODE45 routine
% Loop 1 for Upstream Anechoic p+(0) = 0 sigma(0) = 0
% Initial Conditions
p_plus_in_a = 0; 
p_minus_in_a = 1;   
epsilon_in_a = 0;
ynew0 = [p_plus_in_a p_minus_in_a epsilon_in_a];
% Calculation
x_star_num = xspan;
optsnew = odeset('RelTol',1e-7,'AbsTol',1e-10);
[~,ynew] = ode45(@(xnew,ynew) odex3LEEfcnNonIsentropic_ScatteringMatrix(xnew,ynew,xspan,M,Gamma,Beta_star, ko), x_star_num, ynew0, optsnew);

p_plus_out_a = ynew(end,1); p_minus_out_a = ynew(end,2); epsilon_out_a = ynew(end,3);

% Loop 2 for p+(0) = 0 p-(0) = 0
% Initial Conditions
p_plus_in_b = 0; 
p_minus_in_b = 0;
epsilon_in_b = 1;
ynew0 = [p_plus_in_b p_minus_in_b epsilon_in_b];
% Calculation
optsnew = odeset('RelTol',1e-7,'AbsTol',1e-10);                        
[~,ynew] = ode45(@(xnew,ynew) odex3LEEfcnNonIsentropic_ScatteringMatrix(xnew,ynew,xspan,M,Gamma,Beta_star, ko), x_star_num, ynew0, optsnew);

p_plus_out_b = ynew(end,1); p_minus_out_b = ynew(end,2); epsilon_out_b = ynew(end,3);

% Loop 3 for p-(0) = 0 epsilon(0) = 0
% Initial Conditions
p_plus_in_c = 1; 
p_minus_in_c = 0;
epsilon_in_c = 0;
ynew0 = [p_plus_in_c p_minus_in_c epsilon_in_c];
% Calculation
optsnew = odeset('RelTol',1e-7,'AbsTol',1e-10);
[~,ynew] = ode45(@(xnew,ynew) odex3LEEfcnNonIsentropic_ScatteringMatrix(xnew,ynew,xspan,M,Gamma,Beta_star, ko), x_star_num, ynew0, optsnew);

p_plus_out_c = ynew(end,1); p_minus_out_c = ynew(end,2); epsilon_out_c = ynew(end,3);

% Coefficients of Scattering Matrix :
S_12_3LEE = p_plus_out_a/p_minus_out_a; 
S_21_3LEE = -S_12_3LEE * p_minus_out_c/p_plus_in_c;
S_22_3LEE = p_minus_in_a/p_minus_out_a;

S_13_3LEE = (1/epsilon_out_b) * (p_plus_out_b - S_22_3LEE * p_minus_out_b); 
S_11_3LEE = (p_plus_out_c - S_22_3LEE * p_minus_out_c)/p_plus_in_c;
S_23_3LEE = -S_12_3LEE * p_minus_out_b/epsilon_out_b;

S_32_3LEE = epsilon_out_a/p_minus_out_a;
S_31_3LEE = (epsilon_out_c - S_32_3LEE * p_minus_out_c)/p_plus_in_c;
S_33_3LEE = (1/epsilon_out_b) * (epsilon_out_b - S_32_3LEE * p_minus_out_b);

S_3LEE = [S_11_3LEE S_12_3LEE S_13_3LEE; S_21_3LEE S_22_3LEE S_23_3LEE; S_31_3LEE S_32_3LEE S_33_3LEE];
end

function [M, Beta_star] = MeanFlowForScatteringMatrixNonIsentropicLinear(b_star, M_in, scale) 
% Geometrical description
x_star_num = (0:scale:1)';
 % Temperature Distribution
Beta_star = b_star ./( 1 + b_star .* x_star_num);                    % Normalised Temperature gradient
P_star_in = 1;
% Mean flow calculation
x_span_star = (0:scale:1)'; y0 = [M_in P_star_in];
opts = odeset('RelTol',1e-7,'AbsTol',1e-10);
[~,y] = ode45(@(x,y) odexMeanNonIsentropicfcn(x,y,Beta_star, x_star_num), x_span_star, y0, opts);

% Mean flow Profile
M = y(:,1);
end

function dyprimedx = odex3LEEfcnNonIsentropic_ScatteringMatrix(xnew,ynew,xspan,M,Gamma,Beta_star,He_Omega)
M = interp1(xspan,M,xnew);                           % Interpolate the data set 
He_Omega = interp1(xspan,He_Omega,xnew);
Beta_star = interp1(xspan,Beta_star,xnew);

% P+^ - y(1); P-^ - y(2) sigma^ - y(3)

% dp+/dx = C11*P+ + C12*P- + I_1*sigma^ 
A = ((1i*He_Omega /M) + (Gamma^2)*(M^2)*((Beta_star)/(1 - Gamma*M^2)) + Gamma*Beta_star);
BM_Rhoc = (M^2)*(- Beta_star)/(1 - Gamma*M^2);
I_1 = (1/(2*(1 + M))) * (M^2) * ((Beta_star)/(1-Gamma*(M^2)));

CM_Zeta = - ((Beta_star/2) - (Beta_star)/(1 - Gamma*M^2)); % d/dx(1/Rho c)

% dp-/dx = C21*P+ + C22*P- + I_2*sigma^ 
D = (Beta_star) * (M^2) /(1 - Gamma * (M^2));
E_Rhoc = M * ((1i*He_Omega/M) + (Beta_star)/(1-Gamma*(M^2)));
F_Zeta = -M*((Beta_star/2) - (Beta_star)/(1 - Gamma*M^2));
I_2 = (1/(2*(1 - M))) * (M^2) * ((Beta_star)/(1-Gamma*(M^2)));

% depsilon^ /dx = C31*P+ + C32*P- C33*sigma^ 
Func_M = (- Beta_star*(1 - (M^2)))/(1-Gamma*(M^2));

C_11 = ((-1/(2*(1 + M))) * ( (D + A*M) + (E_Rhoc + BM_Rhoc) + (F_Zeta + CM_Zeta)));
C_12 =  (-1/(2*(1 + M))) * ( (D + A*M) - (E_Rhoc + BM_Rhoc) - (F_Zeta + CM_Zeta));
C_13 = I_1;

C_21 =  (-1/(2*(1 - M))) * ( (D - A*M) + (E_Rhoc - BM_Rhoc) + (F_Zeta - CM_Zeta));
C_22 = ((-1/(2*(1 - M))) * ( (D - A*M) - (E_Rhoc - BM_Rhoc) - (F_Zeta - CM_Zeta)));
C_23 = I_2;

C_31 = Func_M * (Gamma + (1/M));
C_32 = Func_M * (Gamma - (1/M));
C_33 = -(1i*(He_Omega/M) + Gamma*(M^2)*(Beta_star)/(1-Gamma*(M^2)));

dyprimedx = zeros(2,1);
  dyprimedx(1) = C_11* ynew(1) + C_12* ynew(2) + C_13* ynew(3);
  dyprimedx(2) = C_21* ynew(1) + C_22* ynew(2) + C_23* ynew(3);
  dyprimedx(3) = C_31* ynew(1) + C_32* ynew(2) + C_33* ynew(3);
end

% Using ode45 to solve for the mean flow
function dydx = odexMeanNonIsentropicfcn(x,y,Beta_star, x_star_num)
% y(1) - Mach number; y(2) - Pressure (normalised)
% Evaluates the mean flow Mach number and non-dimensional pressure from
% mean conservation first order ODEs
% global Gamma
Gamma = 1.4;                   % Ratio of specific heats for air, treated as a constant

Beta_star = interp1(x_star_num, Beta_star, x);

dydx = zeros(2,1);
dydx(1) = y(1,1) .* ((Beta_star./2) + (Beta_star).*Gamma.*(y(1,1).^2) + (Beta_star).*(Gamma.^2) .* (y(1,1).^4));         
dydx(2) = y(2,1) .* ( - Gamma .* (y(1).^2) .* ((Beta_star).*( 1 + Gamma .* (y(1).^2))));
end