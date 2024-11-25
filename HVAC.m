clc;clear;
%%
variables = {'T', 'T_h', 'T_ext', 'T_env', 'm_dot', 'v', 'c_air', 'A', 'k', 'k_h', 'C', 'C_h', 'u_h', 'h'};

% Define all variables as symbolic using the list
syms(variables{:});

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
dThdt_func = matlabFunction(dThdt, 'Vars', {'T', 'T_h', 'T_env', 'u_h', 'v', 'h', 'k_h', 'm_dot', 'c_air', 'C_h'});

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