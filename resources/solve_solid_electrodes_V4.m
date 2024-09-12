%% SPMe Code (for LG M50 cell)
% Copywrite @Jagmohan
% 2022-02-19

%% This funcion requires 2 inputs and gives 10 output parameters (5 each for n and p electrode)
% Inputs: 2 Parameters
%   1. SoC (State of Charge) and 
%   2. time step of simulation (0.2s for FUDS and 1s for constant and Pulse current)

% Outputs: 10 Parameters
%   1. n-electrode: discrete model(Ad,Bd,Cd,Dd) and initial condition (cse_n_init)
%   2. p-electrode: discrete model(Ad,Bd,Cd,Dd) and initial condition (cse_p_init)

%% Continuous to Discrete time function for solid electrode

function [Ad_cse_n,Bd_cse_n,Cd_cse_n,Dd_cse_n,cse_n_init,Ad_cse_p,Bd_cse_p,Cd_cse_p,Dd_cse_p,cse_p_init] ...
    = solve_solid_electrodes_V4(soc,tstep)
%% Parameters
global D_s_n R_s_n D_s_p R_s_p
% Note: 1. Updated all theta parameters based on OCV parameter results
%       2. Uppdated max concentration values, cse_n_max and cse_p_max
%       3. Updated cse_n_max & cse_p_max vales from LGM50T code
C_bat = 5; % previously 25.67, [Ah] nominal battery capacity (note: my LGM50T cell in the lab has a capacity of 5.0032 Ah)
L_n = 85.2e-6; % [m] thickness of negative electrode
L_sep = 12e-6; % [m] thickness of separator
L_p = 75.6e-6; % [m] thickness of positive electrode
R_s_n = 5.86e-6; % [m] radius of solid negative electrode particle
R_s_p = 5.22e-6; % [m] radius of solid positive electrode particle
Area_n = 0.1027; % [m^2] negative electrode current collector area
Area_p = 0.1027; % [m^2] positive electrode current collector area

%Note: Need to update all epsilon parameters based on the cell
epsilon_e_n = (0.335+0.25)/2; % [-] active material volume fraction (porosity) of electrolyte in negative electrode
epsilon_e_sep = 0.47; % [-] active material volume fraction (porosity) of electrolyte in separator
epsilon_e_p = (0.335+0.25)/2; % [-] active material volume fraction (porosity) of electrolyte in positive electrode
epsilon_s_n = 0.7184; % [-] active material volume fraction of solid negative electrode
epsilon_s_p = 0.5616; % [-] active material volume fraction of solid positive electrode
F = 96485.33289; % [C/mol] Faraday constant
t_plus = 0.2594; % [-] transference number
D_e = 1.7694e-10; % [m^2/s] electrolyte diffusion coefficient (assumed constant for all phases and x-coordinates)
D_s_n = 3.3e-14; % [m^2/s] solid-phase diffusion coefficient for negative electrode
D_s_p = 4e-15; % [m^2/s] solid-phase diffusion coefficient for positive electrode
R = 8.314472; % [J/mol-K] universal gas constant
T = 298.15; % [K] cell temperature
gamma = 0; % [-] activity coefficient
k_n = 6.7160e-12; % [(mol/(s*m^2))*(mol/m^3)^(-3*alpha)] reaction rate in negative electrode
k_p = 3.5446e-11; % [(mol/(s*m^2))*(mol/m^3)^(-3*alpha)] reaction rate in positive electrode
alpha = 0.5; % [-] charge transfer coefficient (for positive and negative electrodes)
R_f_n = 0; % [Ohms*m^2] resistivity of SEI layer for negative electrode
R_f_p = 0; % [Ohms*m^2] resistivity of SEI layer for positive electrode
R_c = 0; %  [Ohms*m^2] contact resistance (resistance of current collectors and wiring)
kappa = 0.9487; % [1/Ohms-m] electrolyte conductivity
% Computed Parameters
a_s_n = 3*epsilon_s_n/R_s_n; % [1/m] negative electrode specific interfacial surface area
a_s_p = 3*epsilon_s_p/R_s_p; % [1/m] positive electrode specific interfacial surface area
epsilon_e = mean([epsilon_e_n,epsilon_e_p]); % [-] uniform active material volume fraction of electrolyte in electrodes (for Pade approximation)
% Limits
V_out_max = 4.2; % [V] terminal voltage upper limit
V_out_min = 2.5; % [V] terminal voltage lower limitSOC_max = 1; % [-] SOC upper limit
SOC_max = 1;
SOC_min = 0; % [-] SOC lower limit
cse_n_max = 33133; % [mol/m^3] max concentration in negative electrode
cse_p_max = 63104; % [mol/m^3] max concentration in positive electrode

%Note: Need to update all theta parameters based on the cell
theta_n_0 = 0.0282; % [-] negative electrode surface stoichiometry when SOC = 0
theta_n_1 = 0.9014; % [-] negative electrode surface stoichiometry when SOC = 1
theta_p_0 = 0.9369; % [-] positive electrode surface stoichiometry when SOC = 0
theta_p_1 = 0.2752; % [-] positive electrode surface stoichiometry when SOC = 1

% Initial Conditions
SOC_0 = soc; % initial SOC
%SOC_0 = 1;
ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte
theta_n(1) = (theta_n_1-theta_n_0)*SOC_0+theta_n_0; % [-] initial negative electrode surface stoichiometry
theta_p(1) = theta_p_0-(theta_p_0-theta_p_1)*SOC_0; % [-] initial positive electrode surface stoichiometry
cse_0_n = theta_n(1)*cse_n_max; % [mol/m^3] initial concentration of ions in negative solid electrode
cse_0_p = theta_p(1)*cse_p_max; % [mol/m^3] initial concentration of ions in positive solid electrode

%% Solving for coefficients: negative electrode
sympref('FloatingPointOutput',true)
format longg

% Note: R_s_n and D_s_n parameters are taken from latest parameters for my
%       cell in lab

% Note: Use below commands to solve for pade function
% Z = solvePade('n');
% c0_n = subs(Z.c0,{'R_s_n','D_s_n'},{5.86e-6,3.3e-14})
% c1_n = subs(Z.c1,{'R_s_n','D_s_n'},{5.86e-6,3.3e-14})
% c2_n = subs(Z.c2,{'R_s_n','D_s_n'},{5.86e-6,3.3e-14})
% d1_n = subs(Z.d1,{'R_s_n','D_s_n'},{5.86e-6,3.3e-14})
% d2_n = subs(Z.d2,{'R_s_n','D_s_n'},{5.86e-6,3.3e-14})

% Solved Pade's Approximation
c0_n = -3/R_s_n;
c1_n = -(4*R_s_n)/(11*D_s_n);
c2_n = -(R_s_n^3)/(165*D_s_n^2);
d1_n = (3*R_s_n^2)/(55*D_s_n);
d2_n = (R_s_n^4)/(3465*D_s_n^2);


% Filling state matrix
A1 = [0,1,0;0,0,1;0,-1/d2_n,-d1_n/d2_n]; 
% Filling the input matrix
B1 = [0;0;1];
% Filling Output matrix
C1 = [c0_n/d2_n,c1_n/d2_n,c2_n/d2_n];
% Filling the Feedthrough Matrix
D1 = 0;

%% n-Electrode: Solving the state space model
% initial concentration for n-electrode
cse_n_init = zeros(3,1); 
cse_n_init(1) = cse_0_n/(c0_n/d2_n);
% time descretization
dt = tstep;

% Discrete time model: n-electrode
% Note: 1. Use below command to use exact solver based on matrix exponent.
%       2. More accurate
% [Ad_cse_n,Bd_cse_n] = c2d(A1,B1,dt);

% Note: Use below command to use Euler integration. Less accurate
Ad_cse_n = A1*dt+eye(size(A1)); % calculated with forward differencing
Bd_cse_n = B1*dt; % calculated with forward differencing

Cd_cse_n = C1 ; Dd_cse_n = D1;

%% Solving for coefficients: positive electrode
sympref('FloatingPointOutput',true)

% Note: R_s_p and D_s_p parameters are taken from latest parameters for my
%       cell in lab

format longg

% Note: Use below commands to solve for pade function
% Z = solvePade('p');
% c0_p = subs(Z.c0,{'R_s_p','D_s_p'},{5.22e-6,4e-15});
% c1_p = subs(Z.c1,{'R_s_p','D_s_p'},{5.22e-6,4e-15});
% c2_p = subs(Z.c2,{'R_s_p','D_s_p'},{5.22e-6,4e-15});
% d1_p = subs(Z.d1,{'R_s_p','D_s_p'},{5.22e-6,4e-15});
% d2_p = subs(Z.d2,{'R_s_p','D_s_p'},{5.22e-6,4e-15});

% Solved Pade's Approximation 
c0_p = -3/R_s_p;
c1_p = -(4*R_s_p)/(11*D_s_p);
c2_p = -(R_s_p^3)/(165*D_s_p^2);
d1_p = (3*R_s_p^2)/(55*D_s_p);
d2_p = (R_s_p^4)/(3465*D_s_p^2);

% Creating the controllable canonical form SS model
% Filling state matrix
A2 = [0,1,0;0,0,1;0,-1/d2_p,-d1_p/d2_p];
% Filling the input matrix
B2 = [0;0;1];
% Filling Output matrix
C2 = [c0_p/d2_p,c1_p/d2_p,c2_p/d2_p];
% Filling the Feedthrough Matrix
D2 = 0;

%% p-Electrode: Solving the state space model
%initial concentration for p-electrode
cse_p_init = zeros(3,1); 
cse_p_init(1) = cse_0_p/(c0_p/d2_p);
% time descretization
dt = tstep;

% Discrete time model
% Note: 1. Use below command to use exact solver based on matrix exponent.
%       2. More accurate
% [Ad_cse_p,Bd_cse_p] = c2d(A2,B2,dt);

% Note: Use below command to use Euler integration. Less accurate
Ad_cse_p = A2*dt+eye(size(A2)); % calculated with forward differencing
Bd_cse_p = B2*dt; % calculated with forward differencing

Cd_cse_p = C2 ; Dd_cse_p = D2;
end




%% Function to get coefficients for Pade's TF
function Z = solvePade(k)

syms c0 c1 c2 d1 d2 m D_s_n R_s_n D_s_p R_s_p 
for i = 1:5
    sprintf('Solving for the %dth derivative',i-1)
    A(i) = mypade_order(i-1);
    B(i) = solveTTF(i-1,k);
end
Z =solve( A(1)==B(1), A(2)==B(2),A(3)==B(3),A(4)==B(4), A(5)==B(5));
end

%% Defining function to calculate Transcendental TF
function [x] = solveTTF(a,k)
%global D_s_n R_s_n
syms TTF(s) Y(s) R_s_n D_s_n R_s_p D_s_p
order = a;
if k=='n'
    Nr(s) = R_s_n*(1-exp(2*R_s_n*sqrt(s/D_s_n))); 
    Dr(s) = D_s_n*(1+R_s_n*sqrt(s/D_s_n)+(R_s_n*sqrt(s/D_s_n)-1)*exp(2*R_s_n*sqrt(s/D_s_n)));
else
    Nr(s) = R_s_p*(1-exp(2*R_s_p*sqrt(s/D_s_p))); 
    Dr(s) = D_s_p*(1+R_s_p*sqrt(s/D_s_p)+(R_s_p*sqrt(s/D_s_p)-1)*exp(2*R_s_p*sqrt(s/D_s_p)));
end
f(s) = s*Nr(s)/Dr(s);
df(s) = diff(f(s),s,order);
[Nr(s),Dr(s)] = numden(df);  % Extracting numerator and denominator from the function
a=0;
while Nr(0)==0 && Dr(0)==0
      dN(s) = diff(Nr,s,1);
      dD(s) = diff(Dr,s,1);
      f(s) = dN(s)/dD(s);
      [Nr(s),Dr(s)] = numden(f(s));
      a= a+1;
      sprintf('running the L''Hospital''s rule for %dth time',a)
end
x = Nr(0)/Dr(0);
end

%% Modified Pade's functions
function P = mypade_order(a)
syms sP(s) c0 c1 c2 d1 d2
sP(s) = (c0 +c1*s+c2*s^2)/(1+d1*s+d2*s^2);
sP1(s) = diff(sP(s),s,a);
P = sP1(0);
end