%% HVAC system university project.
% The protorype is a 20/20/20 cm acrylic box, with some sensors and
% actuators. It is a test environment to test control engineering
% techniques with multiple inputs and outputs.
% The box has an air pipe leading out of the box. On that pipe there is
% first a pair of valves operated by a single motor. This valve allow some
% of the air to be vented to the environment and be replaced by external
% air. If the valve is closed, all air is recirculated, if open, all air is
% vented and external air is drawed in.
% The next component in the pipe is a temperature sensor that measures the
% mixture temperature. This is follower by a heating element, that can
% either be turned on or off, and heats up the mixture. This mixted and
% heated air is then pumped back into the room/box by a motor. This motor's
% speed can be set. Finally, the temperature of this air can also be
% measured.
% There is a temperature sensor inside the box, this is the temperature we
% ultimately want to control. There is also a light sensor in the room
% alongside with an old light bulb. We can also control the light level in
% the room, but we need to pay attention to the fact that the bulb also
% heats up the room.

%% Preparing symbolic toolbox, Parameter values
clc;clear;

% System Constants for Room
% a*b*c sized acrylic box
a = 0.20;                                                                   % side of the box (m)
b = 0.20;                                                                   % side of the box (m)
c = 0.20;                                                                   % side of the box (m)

air_density = 1.225;                                                        % kg/m3 
V = a * b * c;                                                              % m3 air
m_air = V * air_density;                                                    % kg air
                                                    
c_air = 1005;                                                               % Specific heat capacity of air (J/kg·°C)
C_r = m_air * c_air;                                                          % Thermal capacitance of air in the box (J/°C)

% tube calc
d_pipe = 0.035; % m
A_pipe = ((d_pipe/2).^2)*pi;
V_dot_multiplier = A_pipe;
m_dot_multiplier = V_dot_multiplier * air_density;

ss_plant = sys( ...    
    struct( ...  % params
        "T_ext",                20, ...
        "m_dot_multiplier",     m_dot_multiplier, ...
        "Cr",                   C_r, ...
        "Ch",                   3, ...
        "c_air",                c_air, ...
        "k_walls",              2, ...                     % heat loss from room to environment
        "k_heater",             0.5, ...                     % heat loss from heater to environment
        "k_h",                  10, ...                     % heat transfer coefficient from heater to passing air
        "P_l",                  20, ...
        "P_h",                  112.5, ...
        "tau_s",                2, ...
        "tau_v",                2 ...
    ), ...
    {'u_s', 'u_v', 'u_l', 'u_h'}, ...                        % U
    {'T', 'T_h', 's', 'v'}, ...                              % X
    {'T', 'T_hi', 'T_ho', 's', 'v', 'l'} ...                 % Y
);

ss_plant.X_init = [20;20;0.2;0];


clearvars -except ss_plant
%% define equations 
% Define the symbolic equations using the variables from the system object
S = ss_plant.symbols; % Access symbolic variables
s = S.s + 0.2;
m_dot = s * S.m_dot_multiplier;

T_hi = S.v * S.T + (1 - S.v) * S.T_ext;
T_ho = T_hi + S.k_h * S.u_h / s; % Prevent division by zero
l = 1.3 + 1.6 * S.u_l;


dT = ( ...
    S.k_walls * (S.T_ext - S.T) ...
    + S.P_l * S.u_l ...
    + m_dot * S.c_air * (T_ho - S.T) ...
) / S.Cr;

dT_h = ( ...
    S.k_heater * (S.T_ext - S.T_h) ...
    + S.k_h * m_dot * S.c_air * (T_hi - T_ho) ...
    + S.P_h * S.u_h ...
) / S.Ch;

ds = (S.u_s - S.s) / S.tau_s;
dv = (S.u_v - S.v) / S.tau_v;



% defining the system
ss_plant = ss_plant.defineDynamics([dT, dT_h, ds, dv], [S.T, T_hi, T_ho, S.s, S.v, l]);
clearvars -except ss_plant S

%% solving for X, which is STUPID
% f_eqn = sym.empty;
% h_eqn = sym.empty;
% for i=1:length(ss_plant.X)
%     f_eqn(i) = ss_plant.X(:,i) == int(ss_plant.f(:,i), ss_plant.X(:,i));
% end
% for i=1:length(ss_plant.Y)
%     h_eqn(i) = ss_plant.Y(:,i) == ss_plant.h(:,i);
% end
% [f_eqn,h_eqn]
% solve(f_eqn)
% solve([f_eqn,h_eqn], ss_plant.X)
% ss_plant.h
% ss_plant.U
%
%ss_plant.X
%[diff(S.T), diff(S.T_h), diff(S.s), diff(S.v)]
% 
% ss_plant.Y
% solve([[diff(S.T), diff(S.T_h), diff(S.s), diff(S.v)] == ss_plant.f; ...
%     ss_plant.Y == ss_plant.h], ss_plant.X)
% %% compute model properties
% plant = ss_plant.toMatlabFunction("plantFunction");
 

%% linearizing
% [X_dot, Y] = plant([20, 20, 20], [20, 0.01, 0.01, 1, 0]);

ss_plant = ss_plant.linearize( ...
    [S.T, S.T_h, S.s, S.v], ...
    [20, 70, 0.8, 0.5] ...
    );


ss_plant = ss_plant.computeSS();
 mPretty(ss_plant.A);
 mPretty(ss_plant.B);
 mPretty(ss_plant.C);
 mPretty(ss_plant.D);
 model = ss(double(ss_plant.getSS()));
 
 isstable(model)
 step(model)