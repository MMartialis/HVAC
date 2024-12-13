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
a = 0.20;                             % side of the box (m)
b = 0.20;                             % side of the box (m)
c = 0.20;                             % side of the box (m)
rho = 1.225;                          % kg/m3 
V = a * b * c;                        % m3 air
m_air = V * rho;                      % kg air
A_room = 2*a*b + 2*b*c + 2*a*c;       % Area of the walls of the room (m^2)
c_air = 1005;                         % Specific heat capacity of air (J/kg·°C)
C_r = m_air * c_air;                  % Thermal capacitance of air in the box (J/°C)

% tube calc
d_pipe = 0.035;                       % diameter of the pipe (m)
A_pipe = ((d_pipe/2).^2)*pi;          % circumference of the pipe (m) 

ss_plant = sys( ...    
    struct( ...  % params
        "T_ext",                20, ...
        "A_pipe",               A_pipe, ...
        "A_room",               A_room, ...
        "rho",                  rho, ...
        "Cr",                   C_r, ...
        "Ch",                   3, ...
        "c_air",                c_air, ...
        "k_walls",              5.6, ...                     % heat loss from room to environment
        "k_heater",             0.1, ...                     % heat loss from heater to environment
        "P_l",                  20, ...
        "P_h",                  112.5, ...
        "tau_s",                2, ...
        "tau_v",                2, ...
        "tau_l",                0.02, ...
        "s_min",                0.2, ...
        "l_min",                1.3, ...
        "l_max",                8 ...
    ), ...
    {'u_s'; 'u_v'; 'u_l'; 'u_h'}, ...    % U
    {'T'; 'T_h'; 's'; 'v'; 'l'}, ...     % X
    {'T'; 'T_ho'; 's'; 'v'; 'l'} ...     % Y
);

ss_plant.X_init = [21;21;0.2;0.5;1.3];

clearvars -except ss_plant
%% define equations 
% Define the symbolic equations using the variables from the system object
S = ss_plant.symbols; % Access symbolic variables
s = S.s + S.s_min;
m_dot = s * S.rho * S.A_pipe;

T_hi = S.v * S.T + (1 - S.v) * S.T_ext;
T_ho = T_hi + (S.c_air * S.rho * s * S.A_pipe * (S.T_h - T_hi) ...
        / (m_dot * S.c_air));

dT = ( ...
    S.k_walls * S.A_room * (S.T_ext - S.T) ...
    + S.P_l * S.u_l ...
    + m_dot * S.c_air * (T_ho - S.T) ...
) / S.Cr;

dT_h = ( ...
    S.k_heater * (S.T_ext - S.T_h) ...
    - S.c_air * S.rho * s * S.A_pipe * (S.T_h - T_hi) ...
    + S.P_h * S.u_h ...
) / S.Ch;

ds = (S.u_s - S.s) / S.tau_s;
dv = (S.u_v - S.v) / S.tau_v;
dl = ((S.l_min + S.u_l * (S.l_max - S.l_min)) - S.l) / S.tau_l;

% defining the system
ss_plant = ss_plant.defineDynamics( ...
    [dT; dT_h; ds; dv; dl], ...                             % f
    [S.T; T_ho; S.s + S.s_min; S.v; S.l + S.l_min], ...     % h
    [S.T; S.T_h; S.s; S.v; S.l] ...                         % linVars
    );
% Creating function for nonlinear plant block
ss_plant.toMatlabFunction("plantFunction");
% How much faster should the observer be compared to the plant
ss_plant.observerMultiplier = 5;

ss_plant.Q = diag([1 0.001 0.001 0.001 3]);     % {'T'; 'T_h'; 's'; 'v'; 'l'}
ss_plant.R = diag([0.001 0.001 3 1]);           % {'u_s'; 'u_v'; 'u_l'; 'u_h'}
% list of linearization points for creating the mesh
ss_plant.breakpoints = {20:20:80, 20:20:200, 0.2:1:6.2, 0:0.25:1, 1:2:6};
% generating mesh file. This will run for a while :(
ss_plant = ss_plant.getMesh();

