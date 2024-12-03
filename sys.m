classdef sys
    %SYS Class for symbolic system definition
    %   Creates symbolic variables and equations for dynamic systems.
    
    properties
        params      % Struct containing parameter names and values
        symbols     % Struct to store symbolic variables
        X_init      % Initial conditions for plant simulation
        U           % Cell array of input names
        X           % Cell array of state names
        Y           % Cell array of output names
        f           % Array of differential equations (symbolic)
        h           % Output equations (symbolic)
        h_X
        h_U
        h_X_num
        h_U_num
        f_num       % Array of differential equations (symbolic)
        h_num       % Output equations (symbolic)
        f_lin       % Linearized state equations (symbolic)
        h_lin       % Linearized output equations (symbolic)
        f_lin_num   % Linearized state equations with parameters (symbolic)
        h_lin_num   % Linearized output equations with parameters (symbolic)
        A
        B
        C
        D
        A_sym
        B_sym
        C_sym
        D_sym
        A_num
        B_num
        C_num
        D_num
        L
        K
        N
    end

    methods
        function obj = sys(params, U, X, Y)
            % Constructor to initialize parameters, inputs, states, and outputs
            obj.params = params;

            % Create symbolic variables
            paramNames = fieldnames(params);
            symsList = [paramNames; U(:); X(:); Y(:)];
            syms(symsList{:}); % Create symbolic variables

            % Store symbolic variables in the 'params' property
            obj.symbols = struct();
            for i = 1:numel(symsList)
                obj.symbols.(symsList{i}) = eval(symsList{i});
            end
            obj.U = arrayfun(@(x) obj.symbols.(x{1}), U, 'UniformOutput', true);
            obj.X = arrayfun(@(u) obj.symbols.(u{1}), X, 'UniformOutput', true);
            obj.Y = arrayfun(@(y) obj.symbols.(y{1}), Y, 'UniformOutput', true);
            
        end

        function obj = defineDynamics(obj, f, h)
            %DEFINE DYNAMICS Define state and output equations
            %   f: symbolic array of state equations
            %   h: symbolic array of output equations
            
            % Check dimensions
            if numel(f) ~= numel(obj.X)
                error('Number of state equations (f) must match number of states (X).');
            end
            if numel(h) ~= numel(obj.Y)
                error('Number of output equations (h) must match number of outputs (Y).');
            end

            obj.f = f;
            obj.h = h;

            % Substitute parameter values
            paramNames = fieldnames(obj.params);
            paramValues = struct2cell(obj.params);

            obj.f_num = subs(obj.f, paramNames, paramValues);
            obj.h_num = subs(obj.h, paramNames, paramValues);
        end

        function funcHandle = toMatlabFunction(obj, filename)
            %TO MATLAB FUNCTION Generate a function for Simulink or numerical use
            
            % Ensure numerical equations are defined
            if isempty(obj.f_num) || isempty(obj.h_num)
                error('Numerical dynamics are not defined. Call defineDynamics first.');
            end

            % Convert numerical state and output dynamics to MATLAB function
            funcHandle = matlabFunction(obj.f_num, obj.h_num, ...
                'Vars', {obj.X, obj.U}, ...
                'Outputs', {'X_dot', 'Y'});

            func_str = sprintf([
                'function [X_dot, Y]= %s(X, U) \n', ...
                'X = X(:).'';\n', ...
                'U = U(:).'';\n', ...
                'plant = %s;\n', ...
                '[a, b] = plant(X, U);\n', ...
                'X_dot = a(:);\n', ...
                'Y = b(:);\n' ...
                'end'], filename, func2str(funcHandle));

            % Write to file
            fid = fopen(filename + '.m', 'w');
            fprintf(fid, func_str);
            fclose(fid);

        end

        function funcHandle = toObserverFunction(obj, filename)
            % Generates a function handle for nonlinear observer calculation
            % Inputs:
            %   obj: The sys object
            % Outputs:
            %   funcHandle: A function handle that computes X_hat from U and Y
            
            % Ensure symbolic variables exist for X, U, and Y
            if isempty(obj.X) || isempty(obj.U) || isempty(obj.Y)
                error('System must have symbolic X, U, and Y defined.');
            end
            
            % Solve for X in terms of U and Y
            % Equation to solve: Y = h(X, U) -> solve for X
            try
                X_hat_sym = solve(obj.h_num == obj.Y, obj.X);  % Rearrange for X
                % Convert to a flat symbolic array if needed
                if isstruct(X_hat_sym)
                    X_hat_sym = struct2cell(X_hat_sym);  % Convert struct to cell
                end
                X_hat_sym = vertcat(X_hat_sym{:});  % Concatenate into a vector
            catch ME
                error('Could not solve for X. Ensure h is invertible with respect to X.\n%s', ME.message);
            end
        
            % Generate function handle for observer
            funcHandle = matlabFunction(X_hat_sym, ...
                'Vars', {obj.U, obj.Y}, ...
                'Outputs', {'X_hat'});

            func_str = sprintf([
                'function X_hat= %s(U, Y) \n', ...
                'U = U(:).'';\n', ...
                'Y = Y(:).'';\n', ...
                'plant = %s;\n', ...
                'X_hat = plant(X, U);\n', ...
                'end'], filename, func2str(funcHandle));

            % Write to file
            fid = fopen(filename + '.m', 'w');
            fprintf(fid, func_str);
            fclose(fid);
        end


        function obj = linearize(obj, linVars, linVals)
            %LINEARIZE Linearize the system around given points
            %   linPoints: Struct with names and values for linearization

            if isempty(obj.f) || isempty(obj.h)
                error('Define the dynamics (f and h) before linearization.');
            end

            % Use helper function to linearize
            [obj.f_lin, obj.f_lin_num] = obj.performLinearization(obj.f, linVars, linVals);
            [obj.h_lin, obj.h_lin_num] = obj.performLinearization(obj.h, linVars, linVals);
        end
        
        function obj = computeSS(obj)
            %COMPUTESTATESPACE Compute A, B, C, D state-space matrices
            %   Computes symbolic and numerical state-space matrices based on dynamics
            
            if isempty(obj.f) || isempty(obj.h)
                error('Dynamics (f and h) must be defined before computing state-space matrices.');
            end
    
            % Decompose outputs into state- and input-dependent parts
            
            obj.h_X = subs(obj.h_lin, obj.U, zeros(size(obj.U)));                   % Output dependent only on states
            obj.h_U = subs(obj.h_lin, obj.X, zeros(size(obj.X)));                                              % Output dependent only on inputs

            obj.h_X_num = subs(obj.h_lin_num, obj.U, zeros(size(obj.U)));                   % Output dependent only on states
            obj.h_U_num = subs(obj.h_lin_num, obj.X, zeros(size(obj.X)));  

            % Calculate Jacobians for the A, B, C, and D matrices
            obj.A_sym = jacobian(obj.f_lin, obj.X);             % System matrix
            obj.B_sym = jacobian(obj.f_lin, obj.U);             % Input matrix
            obj.C_sym = jacobian(obj.h_X, obj.X);               % Output matrix for state response
            obj.D_sym = jacobian(obj.h_U, obj.U);               % Direct input-output relationship
    
            % Calculate Jacobians for the A, B, C, and D matrices
            obj.A = jacobian(obj.f_lin_num, obj.X);             % System matrix
            obj.B = jacobian(obj.f_lin_num, obj.U);             % Input matrix
            obj.C = jacobian(obj.h_X_num, obj.X);               % Output matrix for state response
            obj.D = jacobian(obj.h_U_num, obj.U);               % Direct input-output relationship
            
            obj.A_num = double(obj.A);
            obj.B_num = double(obj.B);
            obj.C_num = double(obj.C);
            obj.D_num = double(obj.D);

        end
    
        function [A, B, C, D] = getSS(obj)
            %GETSTATESPACE Return numerical state-space matrices
            %   Outputs:
            %       A, B, C, D: Numerical state-space matrices
    
            if isempty(obj.A)
                error('State-space matrices have not been computed yet. Call computeSS first.');
            end
    
            A = obj.A;
            B = obj.B;
            C = obj.C;
            D = obj.D;
        end
    end

    methods (Access = private)
        function [linEq, linNum] = performLinearization(obj, eq, subsVars, subsVals)
            %PERFORMLINEARIZATION Helper to linearize equations
            %   eq: Array of symbolic equations to linearize
            %   subsVars: Cell array of variables to linearize
            %   subsVals: Struct with values for linearization

            % Substitute the linearization points

            % Linearize using Taylor expansion
            linEq = taylor(eq, subsVars, 'ExpansionPoint', subsVals, 'Order', 2);

            % Substitute parameter values into the linearized model
            paramNames = fieldnames(obj.params);
            paramValues = struct2cell(obj.params);

            linNum = subs(linEq, paramNames, paramValues);
        end

    end
end
