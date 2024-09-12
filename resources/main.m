
%% SPMe Code (for LG M50 cell)
% Copywrite @Jagmohan
% 2022-02-19

clc;
clear all; close all;

tic  % Starting clock to check for run time
%% Parameters
C_bat = 5.032; % [Ah] nominal battery capacity (note: my LGM50T cell in the lab has a capacity of 5.0032 Ah)
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

epsilon_s_n = 0.7374; % [-] active material volume fraction of solid negative electrode
epsilon_s_p = 0.5758; % [-] active material volume fraction of solid positive electrode

t_plus = 0.2594; % [-] transference number
D_e = 1.7694e-10; % [m^2/s] electrolyte diffusion coefficient (assumed constant for all phases and x-coordinates)
D_s_n = 3.3e-14; % [m^2/s] solid-phase diffusion coefficient for negative electrode
D_s_p = 4e-15; % [m^2/s] solid-phase diffusion coefficient for positive electrode

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
V_out_max = 4.5; % [V] terminal voltage upper limit
V_out_min = 2.5; % [V] terminal voltage lower limit
SOC_max = 1; % [-] SOC upper limit
SOC_min = 0; % [-] SOC lower limit
cse_n_max = 33133; % [mol/m^3] max concentration in negative electrode
cse_p_max = 63104; % [mol/m^3] max concentration in positive electrode

%Note: Need to update all theta parameters based on the cell
theta_n_0 = 0.0282; % [-] negative electrode surface stoichiometry when SOC = 0
theta_n_1 = 0.9014; % [-] negative electrode surface stoichiometry when SOC = 1
theta_p_0 = 0.9369; % [-] positive electrode surface stoichiometry when SOC = 0
theta_p_1 = 0.2752; % [-] positive electrode surface stoichiometry when SOC = 1

R = 8.314472; % [J/mol-K] universal gas constant
T = 298.15; % [K] cell temperature
F = 96485.33289; % [C/mol] Faraday constant

%% Input data
% Note: Change the file name to 'SPMe_results_pulse1C_v1.mat' or 'SPMe_results_FUDS_v1.mat'
%       to simulate other current profiles

% x = load('Pulse_CurrentProfile_5A.mat');  % Importing current inputs 
% x = load('FUDS_CurrentProfile.mat');
x = load('1CC_CurrentProfile_5A.mat');

t = x.time; % simulation time
dt = t(2)-t(1);
I = -x.I; %[A] deintercalation current  (+ve is discharging)
j_n = +I./(Area_n*L_n);   % curren flux at n-electrode: j_n(Li) = I/A*thickness 
j_p = -I./(Area_p*L_p);   % curren flux at p-electrode: j_p(Li) = -I/A*thickness 
inp_cse_n = j_n./(F*a_s_n);    % Input for negative electrode
inp_cse_p = j_p./(F*a_s_p);    % Input for positive electrode
inp_ce_n = I;  % Input for electrolyte at n-electrode
inp_ce_p = I;  % Input for electrolyte at p-electrode

%% Calling solid electrode and electrolyte functions
% Initial Conditions
%SOC_0 = x.SOC_0; % initial SOC 
SOC_0 = 1; % [fraction] initial SOC
ce_0 = 1000; % [mol/m^3] initial concentration of Li-ions in electrolyte
theta_n(1) = (theta_n_1-theta_n_0)*SOC_0+theta_n_0; % [-] initial negative electrode surface stoichiometry
theta_p(1) = theta_p_0-(theta_p_0-theta_p_1)*SOC_0; % [-] initial positive electrode surface stoichiometry
cse_0_n = theta_n(1)*cse_n_max; % [mol/m^3] initial concentration of ions in negative solid electrode
cse_0_p = theta_p(1)*cse_p_max; % [mol/m^3] initial concentration of ions in positive solid electrode

% Solving for SPMe Solid Li Electrode sub-model
[Ad_cse_n,Bd_cse_n,Cd_cse_n,Dd_cse_n,cse_n_init,Ad_cse_p,Bd_cse_p,Cd_cse_p,Dd_cse_p,cse_p_init] ...
    = solve_solid_electrodes_V4(SOC_0,dt);
% Solving for SPMe Electrolyte sub-model
[Ad_ce_n,Bd_ce_n,Cd_ce_n,Dd_ce_n,ce_n_init,Ad_ce_p,Bd_ce_p,Cd_ce_p,Dd_ce_p,ce_p_init] ...
    = solve_electrolyte_V4(SOC_0,dt);

%% Simulation Setup
k = length(t);
n_se = 3; % order of 3rd order Pade's TF in solid electrode 
n_e = 1; % order of 1st order Pade's TF in electrolyte 
% solid electrode
Xn_cse = zeros(n_se,k+1);
Xp_cse = zeros(n_se,k+1);
cse_n = zeros(1,k);
Xn_cse(:,1) = cse_n_init; % initial conc at n-electrode surface
cse_p = zeros(1,k);
Xp_cse(:,1) = cse_p_init; % initial conc at p-electrode surface
% electrolyte
Xn_ce = zeros(n_e,k+1);
Xp_ce = zeros(n_e,k+1);
ce_n = zeros(1,k);
Xn_ce(:,1) = ce_n_init; % initial conc at n-electrolyte
ce_p = zeros(1,k);
Xp_ce(:,1) = ce_p_init; % initial conc at p-electrolyte
% Li intercalation/deintercalation 
d_n = L_n; % [m] n-electrode thickness
d_p = L_p; % [m] p-electrode thickness
d_sep = L_sep; % [m] separator thickness
% effective resistance
R_omega = R_c/Area_n + R_f_n/(Area_n*d_n*a_s_n) + R_f_p/(Area_p*d_p*a_s_p) + ...
    (d_n + d_p)/(2*Area_n*kappa*epsilon_e^1.5) + (d_sep)/(Area_n*kappa*epsilon_e_sep^1.5);
% Li Intercalation/deintercalation        
i_on = zeros(1,k);
i_op = zeros(1,k);
zeta_n = zeros(1,k);
zeta_p = zeros(1,k);
del_phi_e_con = zeros(1,k);
% Open Circuit Potential
theta_n = zeros(1,k);
theta_p = zeros(1,k);
U_n = zeros(1,k);  % n-electrode OCP
U_p = zeros(1,k);  % p-electrode OCP
% Battery terminal voltage and overpotential
eta_n = zeros(1,k);
eta_p = zeros(1,k);
V_out = zeros(1,k);
SOC_n = zeros(1,k); % [-] SOC from negative electrode surface stoichiometry
SOC_p =  zeros(1,k);% [-] SOC from positive electrode surface stoichiometry
SOC = zeros(1,k); % mean SOC

%% Simulating the overall system
for i =1:k  
    
    % solid n-electrode
    Xn_cse(:,i+1) = Ad_cse_n*Xn_cse(:,i) + Bd_cse_n*inp_cse_n(i);
    cse_n(:,i) = Cd_cse_n*Xn_cse(:,i) + Dd_cse_n*inp_cse_n(i);
    
    % solid p-electrode
    Xp_cse(:,i+1) = Ad_cse_p*Xp_cse(:,i) + Bd_cse_p*inp_cse_p(i);
    cse_p(:,i) = Cd_cse_p*Xp_cse(:,i) + Dd_cse_p*inp_cse_p(i);
    
    % n-electrolyte
    Xn_ce(:,i+1) = Ad_ce_n*Xn_ce(:,i) + Bd_ce_n*inp_ce_n(i);
    ce_n(:,i) = Cd_ce_n*Xn_ce(:,i) + Dd_ce_n*inp_ce_n(i) + ce_0; % ce_0 added to compansate for zero IC
    
    % p-electrolyte
    Xp_ce(:,i+1) = Ad_ce_p*Xp_ce(:,i) + Bd_ce_p*inp_ce_p(i);
    ce_p(:,i) = Cd_ce_p*Xp_ce(:,i) + Dd_ce_p*inp_ce_p(i) + ce_0;  % ce_0 added to compansate for zero IC 
    
    % Li intercalation/deintercalation  
    i_on(i) = F*k_n*(ce_0^alpha)*((cse_n_max - cse_n(i))^alpha)*(cse_n(i)^alpha); % corrected with ce_0
    i_op(i) = F*k_p*(ce_0^alpha)*((cse_p_max - cse_p(i))^alpha)*(cse_p(i)^alpha); % corrected with ce_0
    zeta_n(i) = j_n(i)/(2*a_s_n*i_on(i));
    zeta_p(i) = j_p(i)/(2*a_s_p*i_op(i));
    del_phi_e_con(i) = (2*R*T*(1-t_plus)/F)*(1+gamma)*log(ce_p(i)/ce_n(i));   
    
    % Battery terminal voltage and overpotential
    eta_n(i) = (R*T/(alpha*F))*log(zeta_n(i) + sqrt(zeta_n(i)^2 + 1));
    eta_p(i) = (R*T/(alpha*F))*log(zeta_p(i) + sqrt(zeta_p(i)^2 + 1));
   
    % Open Circuit Potential
    theta_n(i) = cse_n(i)/cse_n_max; % [-] negative electrode surface stoichiometry
    U_n(i) = 1.9793*exp(-39.3631*theta_n(i))+0.2482...
            -0.0909*tanh(29.8538*(theta_n(i)-0.1234))...
            -0.04478*tanh(14.9159*(theta_n(i)-0.2769))...
            -0.0205*tanh(30.4444*(theta_n(i)-0.6103)); % [V] open circuit potential for negative electrode
   
    theta_p(i) = cse_p(i) / cse_p_max; % [-] positive electrode surface stoichiometry
    U_p(i) = -0.809*theta_p(i)+4.4875-0.0428*tanh(18.5138*(theta_p(i)-0.5542))-...
            17.7326*tanh(15.789*(theta_p(i)-0.3117))+...
            17.5842*tanh(15.9308*(theta_p(i)-0.312)); % [V] open circuit potential for positive electrode
 
    % Voltage output
    V_out(i) = U_p(i) - U_n(i) + del_phi_e_con(i) + eta_p(i) - eta_n(i) - I(i)*R_omega;  
    
    % SOC
    SOC_n(i) = (theta_n(i) - theta_n_0)/(theta_n_1 - theta_n_0); % [-] SOC from negative electrode surface stoichiometry
    SOC_p(i) = (theta_p(i) - theta_p_0)/(theta_p_1 - theta_p_0); % [-] SOC from positive electrode surface stoichiometry
    SOC(i) = (SOC_n(i) + SOC_p(i))/2; % [-] mean SOC
    
    if V_out(i) < V_out_min || V_out(i) > V_out_max %|| SOC(i) < SOC_min || SOC(i) > SOC_max % if voltage or SOC limit is reached
        break % stop simulation
    end
end

%% Plotting the results

figure(1)
plot(t,I)
xlabel('time [s]')
ylabel('Current Profile [A]')
title('Current Profile over time')

figure(1)
subplot(3,3,1);
plot(t,V_out)
xlabel('time [s]')
ylabel('Output Voltage [V]')
title('Voltage output over time')

subplot(3,3,2);
plot(t,cse_n)
xlabel('time [s]')
ylabel('n-electrode conc [mol/m^3]')
title('Surface Conc at negative electrode over time')

subplot(3,3,3);
plot(t,cse_p)
xlabel('time [s]')
ylabel('p-electrode conc [mol/m^3]')
title('Surface Conc at positive electrode over time')

subplot(3,3,4);
plot(t,ce_n)
xlabel('time [s]')
ylabel('n-electrolyte conc [mol/m^3]')
title('Electrolyte conc at x=0n over time')

subplot(3,3,5);
plot(t,ce_p)
xlabel('time [s]')
ylabel('p-electrolyte conc [mol/m^3]')
title('Electrolyte conc at x=0p over time')

subplot(3,3,6);
plot(t,eta_n)
xlabel('time [s]')
ylabel('Overpotential at n-electrode [V]')
title('Overpotential at n-electode over time')

subplot(3,3,7);
plot(t,eta_p)
xlabel('time [s]')
ylabel('Overpotential at p-electrode [V]')
title('Overpotential at p-electode over time')

subplot(3,3,8);
plot(t,U_n)
xlabel('time [s]')
ylabel('OCP at n-electrode [V]')
title('OCP at n-electode over time')

subplot(3,3,9);
plot(t,U_p)
xlabel('time [s]')
ylabel('OCP at p-electrode [V]')
title('OCP at p-electode over time')

toc    % Stopping clock to check for run time