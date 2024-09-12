%% SPMe Code (for LG M50 cell)
% Copywrite @Jagmohan
% 2022-02-19

%% This funcion requires 1 inputs and gives 10 output parameters (5 each for n and p electrode)
% Inputs: 2 Parameters
%   1. Time step of simulation (0.2s for FUDS and 1s for constant and Pulse current)

% Outputs: 10 Parameters
%   1. n-electrolyte: discrete model(Ad,Bd,Cd,Dd) and initial condition (ce_n_init)
%   2. p-electrolyte: discrete model(Ad,Bd,Cd,Dd) and initial condition (ce_p_init)

%% Continuous to Discrete time function for solid electrode

function [Ad_ce_n,Bd_ce_n,Cd_ce_n,Dd_ce_n,ce_n_init,Ad_ce_p,Bd_ce_p,Cd_ce_p,Dd_ce_p,ce_p_init] ...
    = solve_electrolyte_V4(soc,tstep)

%% Model Parameters
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
%ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte
theta_n(1) = (theta_n_1-theta_n_0)*SOC_0+theta_n_0; % [-] initial negative electrode surface stoichiometry
theta_p(1) = theta_p_0-(theta_p_0-theta_p_1)*SOC_0; % [-] initial positive electrode surface stoichiometry
cse_0_n = theta_n(1)*cse_n_max; % [mol/m^3] initial concentration of ions in negative solid electrode
cse_0_p = theta_p(1)*cse_p_max; % [mol/m^3] initial concentration of ions in positive solid electrode

%% Coefficient of Pade's TF
% Note: Need to update this based on the new L_n, L_sep,and L_p values
sympref('FloatingPointOutput',true)
syms De Es Ee F t0 A 
g1_On = (1.1294e-21*(8.2785e+32*Ee*Es^1.5000 + 9.4081e+30*Ee^1.5000*Es - 1.1854e+32*Ee^2.5000*t0 - 6.6798e+31*Es^2.5000*t0 + 1.1854e+32*Ee^2.5000 + 6.6798e+31*Es^2.5000 - 8.2785e+32*Ee*Es^1.5000*t0 - 9.4081e+30*Ee^1.5000*Es*t0))/(1.7709e+15*A*De*Ee^1.5000*Es^2.5000*F + 2.3730e+16*A*De*Ee^2.5000*Es^1.5000*F);
h1_On = (3.1612e+81*Ee^2.5000*Es^2.5000 + 1.1157e+82*Ee*Es^4 + 3.1666e+80*Ee^4*Es + 1.2796e+80*Ee^1.5000*Es^3.5000 + 2.0231e+82*Ee^3.5000*Es^1.5000 + 1.7668e+81*Ee^5 + 4.6751e+80*Es^5 + 7.0100e+82*Ee^2*Es^3 + 1.7503e+79*Ee^3*Es^2)/(1.5457e+89*De*Ee^(1/2)*Es^5 + 5.6603e+89*De*Ee^3*Es^2.5000 + 2.5670e+91*De*Ee^2.5000*Es^3 + 2.1770e+88*De*Ee^2*Es^3.5000 + 3.6757e+90*De*Ee^4*Es^1.5000 + 3.9869e+90*De*Ee^1.5000*Es^4);
g1_Lc = -(1.1294e-21*(8.6147e+32*Ee*Es^1.5000 + 9.4081e+30*Ee^1.5000*Es - 1.3360e+32*Ee^2.5000*t0 - 5.9271e+31*Es^2.5000*t0 + 1.3360e+32*Ee^2.5000 + 5.9271e+31*Es^2.5000 - 8.6147e+32*Ee*Es^1.5000*t0 - 9.4081e+30*Ee^1.5000*Es*t0))/(1.7709e+15*A*De*Ee^1.5000*Es^2.5000*F + 2.3730e+16*A*De*Ee^2.5000*Es^1.5000*F);
h1_Lc = (3.1246e+81*Ee^2.5000*Es^2.5000 + 1.0005e+82*Ee*Es^4 + 3.3906e+80*Ee^4*Es + 1.0742e+80*Ee^1.5000*Es^3.5000 + 2.1734e+82*Ee^3.5000*Es^1.5000 + 1.9912e+81*Ee^5 + 3.2662e+80*Es^5 + 7.0885e+82*Ee^2*Es^3 + 1.7503e+79*Ee^3*Es^2)/(1.3715e+89*De*Ee^(1/2)*Es^5 + 6.0086e+89*De*Ee^3*Es^2.5000 + 4.1425e+90*De*Ee^4*Es^1.5000 + 2.6712e+91*De*Ee^2.5000*Es^3 + 2.1770e+88*De*Ee^2*Es^3.5000 + 3.8313e+90*De*Ee^1.5000*Es^4);

format shortG

g1_On = double(subs(g1_On,{Ee,Es,A,De,F,t0},{epsilon_e_n,epsilon_e_sep,Area_n,D_e,96485.33289,t_plus}));
h1_On = double(subs(h1_On,{Ee,Es,A,De,F,t0},{epsilon_e_n,epsilon_e_sep,Area_n,D_e,96485.33289,t_plus}));
g1_Lc = double(subs(g1_Lc,{Ee,Es,A,De,F,t0},{epsilon_e_n,epsilon_e_sep,Area_n,D_e,96485.33289,t_plus}));
h1_Lc = double(subs(h1_Lc,{Ee,Es,A,De,F,t0},{epsilon_e_n,epsilon_e_sep,Area_n,D_e,96485.33289,t_plus}));

%% Creating the controllable canonical form SS model
% Initial Condition
ce_n_init = 0;   % IC for electrolyte at n-electrode
ce_p_init = 0;   % IC for electrolyte at p-electrode
dt = tstep;  %1s for Constant and Pulse while 0.2s for FUDS
%dt = 1;

%% Electrolyte model at n-Electrode
A1 = -1/h1_On;
B1 = 1;
C1 = g1_On/h1_On;
D1 = 0;

% Discrete time model
% Note: 1. Use below command to use exact solver based on matrix exponent.
%       2. More accurate
% [Ad_ce_n,Bd_ce_n] = c2d(A1,B1,dt);

% Note: Use below command to use Euler integration. Less accurate
Ad_ce_n = A1*dt+eye(size(A1)); % calculated with forward differencing
Bd_ce_n = B1*dt; % calculated with forward differencing

Cd_ce_n = C1 ; Dd_ce_n = D1;

%% Electrolyte model at p-Electrode
A2 = -1/h1_Lc;
B2 = 1;
C2 = -g1_Lc/h1_Lc;
D2 = 0;

% Discrete time model
% Note: 1. Use below command to use exact solver based on matrix exponent.
%       2. More accurate

% [Ad_ce_p,Bd_ce_p] = c2d(A2,B2,dt);

% Note: Use below command to use Euler integration. Less accurate
Ad_ce_p = A2*dt+eye(size(A2)); % calculated with forward differencing
Bd_ce_p = B2*dt; % calculated with forward differencing

Cd_ce_p = C2 ; Dd_ce_p = D2;
end