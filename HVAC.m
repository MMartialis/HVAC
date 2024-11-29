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
area = 2*a*b+2*a*c+2*b*c;                                                   % Surface area of walls (m²)
                                                    
c_air = 1005;                                                               % Specific heat capacity of air (J/kg·°C)
C = m_air * c_air;                                                          % Thermal capacitance of air in the box (J/°C)
                                                         

sys_p = sys( ...    
    struct( ...  % params
        'A_room',           area, ...                   
        'alpha_T_ext',      0.001, ...                  
        'C',                C, ...                      
        'c_air',            c_air, ...                      % Specific heat of air (J/(kg*K))
        'C_h',              50, ...                         % Thermal capacitance of heating element (J/K)
        'k_room',           1, ...                          % Heat transfer coefficient of walls (W/kgK)
        'k_h',              1, ...                          % Heat transfer coefficient with environment (W/kgK)
        'K_h',              1, ...                          % Heat transfer coefficient between heater and air (W/K)
        'P_h',              100, ...                        % Heating power when on (W)
        'P_bulb',           20, ...                         % Maximum power of bulb (W)
        'eff_bulb',         0.1 ... 
    ), ...  
    {'T_ext','m_dot_i', 'm_dot_e', 'u_bulb', 'u_h'}, ...    % U
    {'T', 'T_h', 'T_ext'}, ...                              % X
    {'T', 'T_vo', 'T_ho', 'Light'} ...                      % Y
);

sys_c = sys( ...    
    struct( ...  % params
        'A_room',           area, ...                   
        'alpha_T_ext',      0.001, ...                  
        'C',                C, ...                      
        'c_air',            c_air, ...                  % Specific heat of air (J/(kg*K))
        'C_h',              50, ...                     % Thermal capacitance of heating element (J/K)
        'k_room',           1, ...                      % Heat transfer coefficient of walls (W/kgK)
        'k_h',              1, ...                      % Heat transfer coefficient with environment (W/kgK)
        'K_h',              1, ...                      % Heat transfer coefficient between heater and air (W/K)
        'P_h',              100, ...                    % Heating power when on (W)
        'P_bulb',           20, ...                     % Maximum power of bulb (W)
        'eff_bulb',         0.1 ...
    ), ...
    {'m_dot_i', 'm_dot_e', 'u_bulb', 'u_h'}, ...        % U
    {'T', 'T_h', 'T_ext'}, ...                          % X
    {'T', 'T_vo', 'T_ho', 'Light'} ...                  % Y
);

clearvars -except sys_p sys_c
%% define equations 
% Define the symbolic equations using the variables from the system object
s = sys_p.symbols; % Access symbolic variables


% output derivatives
% temperature after the valve
m_dot = s.m_dot_e + s.m_dot_i;

T_vo = ((s.m_dot_e * s.T_ext * s.c_air) + (s.m_dot_i * s.T * s.c_air)) ...
     / (m_dot * s.c_air);
% Temperature after the heater
T_ho = (m_dot * s.c_air * T_vo + s.K_h * s.T_h) ...
     / (m_dot * s.c_air + s.K_h); 
% Light level inside the room
Light = s.P_bulb * s.u_bulb * s.eff_bulb;                                               % Light output

%state derivatives
% room
Q_loss = s.k_room * s.A_room * (s.T - s.T_ext);                                  % Heat loss through walls
Q_heat = s.K_h * m_dot * s.c_air * (s.T_h - T_vo) ...
       / (m_dot * s.c_air + s.K_h); % Heater heat
Q_bulb = s.P_bulb * s.u_bulb * (1 - s.eff_bulb);                            % Bulb heat

% heater
Qh_heater = s.u_h * s.P_h;              
Qh_wall = s.k_h * (s.T_h - s.T_ext);                
Qh_vent = m_dot * s.c_air * (s.T_h - T_vo);

% state derivatives
dTdt = (-Q_loss + Q_heat + Q_bulb) / s.C;   
dThdt = (Qh_heater - Qh_wall - Qh_vent) / s.C_h;                                     % Heater temperature rate
dTextdt = s.alpha_T_ext * (s.T - s.T_ext);                                           % External temp rate

% defining the system

sys_p = sys_p.defineDynamics([dTdt, dThdt, dTextdt], [s.T, T_vo, T_ho, Light]);
sys_c = sys_c.defineDynamics([dTdt, dThdt, dTextdt], [s.T, T_vo, T_ho, Light]);
clearvars -except sys_p sys_c s

%% compute model properties
% plant
plant = sys_p.toMatlabFunction();
X_init = [20;20;20];

[X_dot, Y] = plant([20, 20, 20], [20, 0.01, 0.01, 1, 0]);

% controller
sys_c = sys_c.linearize( ...
    [s.T, s.T_h, s.T_ext, s.m_dot_i, s.m_dot_e], ...
    [20, 120, 20, 0.01, 0.01] ...
    );



sys_c = sys_c.computeSS();

 mPretty(sys_c.A);
 mPretty(sys_c.B);
 mPretty(sys_c.C);
 mPretty(sys_c.D);
 model = ss(double(sys_c.getSS()));
 
 isstable(model)
 step(model)