clc;clear;
variables = {'T', 'T_h', 'T_ext', 'T_env', 'm_dot', 'v', 'c_air', 'A', 'k', 'k_h', 'C', 'C_h', 'u_h', 'h'};


% Define all variables as symbolic using the list
syms(variables{:});

% Define state vector and input vector
ss_X = [T; T_h];
ss_U = [u_h; v; m_dot; T_ext];

% Equations for heat change inside the room
Q_loss = -k * A * (T - T_ext);              % heat loss through the walls

Q_vent = m_dot * c_air * v * (T_ext - T);   % heat change caused by the
                                            % ventillated air mixture

Q_heat = m_dot * c_air * v * (T_h - T);     % added heat by heater element

% Symbolic differential equations

% Room temperature change (dT/dt)
dTdt = (Q_loss + Q_vent + Q_heat) / C;
pretty(dTdt);

% Heating element temperature change (dT_h/dt)
dThdt = (u_h * h - k_h * (T_h - T_env) - m_dot * c_air * v * (T_h - T)) / C_h;
pretty(dThdt);

% Define f (state derivatives)
ss_f = [dTdt; dThdt];

clear(variables{:});
%% 
% System Constants for Room
A = 10;           % Surface area of walls (m²)
k = 0.5;          % Heat transfer coefficient of walls (W/m²°C)
C = 1000;         % Thermal capacitance of air in the box (J/°C)
c_air = 1005;     % Specific heat capacity of air (J/kg·°C)
T_ext = 20;        % External temperature (°C)

% Heating Element Constants
C_h = 50;         % Thermal capacitance of heating element (J/°C)
k_h = 0.1;        % Heat transfer coefficient with environment (W/°C)
h = 200;          % Heating power when on (W)
T_env = 20;       % Environmental temperature (°C)

% Control Variables
v = 1;          % Valve position (fraction of external air, 0≤v≤1)
u_h = 0;          % Heating element control (1 for on, 0 for off)
m_dot = 0.1;     % Mass flow rate of air (kg/s)

% Initial Conditions
T_init = 30;      % Initial room temperature (°C)
T_h_init = 30;   % Initial heating element temperature (°C)
initial_conditions = [T_init; T_h_init];

%% Parameter values (stored in a struct for easy reuse)
params = struct( ...
    'A', 10, ...
    'k', 0.5, ...
    'C', 1000, ...
    'c_air', 1005, ...
    'C_h', 50, ...
    'k_h', 0.1, ...
    'h', 200 ...
    );


% % Verify equations are now numeric
% dTdt_numeric = eval(dTdt);
% dThdt_numeric = eval(dThdt);
% 
% disp('Initial Rate of Room Temperature change:');
% disp(dTdt_numeric);
% disp('Initial Rate of Heating Element Temperature change:');
% disp(dThdt_numeric);
% Convert symbolic equations to numeric function handles
dTdt_func = matlabFunction(dTdt, 'Vars', {'T', 'T_h', 'T_ext', 'm_dot', 'v', 'c_air', 'A', 'k', 'C'});
dThdt_func = matlabFunction(dThdt, 'Vars', {'T', 'T_h', 'T_env', 'u_h','v', 'h', 'k_h', 'm_dot', 'c_air', 'C_h'});

% Define the system of equations for ode45 using these functions
dTdt_system = @(t, X) [
    dTdt_func(X(1), X(2), T_ext, m_dot, v, c_air, A, k, C);
    dThdt_func(X(1), X(2), T_env, u_h, v, h, k_h, m_dot, c_air, C_h);
];

% Define time span for simulation
tspan = [0 1000];

% Use ode45 to simulate the temperature evolution over time
[t, X] = ode45(dTdt_system, tspan, initial_conditions);

% Plot results
figure(1);
clf;
plot(t, X(:, 1), 'DisplayName', 'Room Temperature (T)');
hold on;
plot(t, X(:, 2), 'DisplayName', 'Heating Element Temperature (T_h)');
xlabel('Time (s)');
ylabel('Temperature (°C)');
legend;
title('Temperature Evolution in Room and Heating Element');

% Calculate Jacobians symbolically
A_sym = jacobian(ss_f, ss_X);   % System matrix
B_sym = jacobian(ss_f, ss_U);   % Input matrix
C_sym = eye(2);           % Output matrix (identity)
D_sym = zeros(2, 3);      % Direct transmission matrix (zeros)

% Substitute values for symbolic matrices
A_num = subs(A_sym, {"A","k","C","c_air","C_h","k_h","h"}, {params.A,params.k,params.C,params.c_air,params.C_h,params.k_h,params.h});
B_num = subs(B_sym, {"A","k","C","c_air","C_h","k_h","h"}, {params.A,params.k,params.C,params.c_air,params.C_h,params.k_h,params.h});
C_num = double(C_sym);  % Identity matrix
D_num = double(D_sym);  % Zero matrix


% Display symbolic state-space matrices
disp('Symbolic A matrix:');
pretty(A_sym)
disp('Symbolic B matrix:');
pretty(B_sym)
disp('Symbolic C matrix:');
disp(C_sym);
disp('Symbolic D matrix:');
disp(D_sym);


% Display numerical state-space matrices
disp('Numeric A matrix:');
pretty(A_num);
disp('Numeric B matrix:');
pretty(B_num);
disp('Numeric C matrix:');
disp(C_num);
disp('Numeric D matrix:');
disp(D_num);