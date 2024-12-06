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
        Q
        R
        f           % Array of differential equations (symbolic)
        h           % Output equations (symbolic)
        f_num       % Array of differential equations (symbolic)
        h_num       % Output equations (symbolic)
        linVars
        mesh
        breakpoints
        meshSizes
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

        function obj = defineDynamics(obj, f, h, linVars)
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
            obj.linVars = linVars;
            
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
                'X = X(:);\n', ...
                'U = U(:);\n', ...
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

        function SS = getSS(obj, linVals)
            %LINEARIZE Linearize the system around given points
            %   linPoints: Struct with names and values for linearization

            if isempty(obj.f) || isempty(obj.h)
                error('Define the dynamics (f and h) before linearization.');
            end
            % Use helper function to linearize
            [f_lin, f_lin_num] = obj.performLin(obj.f, linVals);
            [h_lin, h_lin_num] = obj.performLin(obj.h, linVals);
            
            % devide the ouptut equition into two parts, those dependent on
            % inputs and those on states
            h_X_num = subs(h_lin_num, obj.U, zeros(size(obj.U)));                   % Output dependent only on states
            h_U_num = subs(h_lin_num, obj.X, zeros(size(obj.X)));  

            % Calculate Jacobians for the A, B, C, and D matrices
            SS.A = double(jacobian(f_lin_num, obj.X));             % System matrix
            SS.B = double(jacobian(f_lin_num, obj.U));             % Input matrix
            SS.C = double(jacobian(h_X_num, obj.X));               % Output matrix for state response
            SS.D = double(jacobian(h_U_num, obj.U));               % Direct input-output relationship
        end

        function SS = getL(obj, SS, multiplier)
            
            Obs_matrix = obsv(SS.A, SS.C);
            if rank(Obs_matrix) < size(SS.A, 1)
                error('System is not fully observable.');
            end

            desired_observer_poles = eig(SS.A) * multiplier;  
            
            SS.L = place(SS.A', SS.C', desired_observer_poles)';  
        end

        function SS = getX(obj, SS, Q, R)
            SS.K = lqr(SS.A, SS.B, Q, R);
        end


        function obj = getMesh(obj)
            % Determine the size of the mesh
            obj.meshSizes = cellfun(@length, obj.breakpoints);  % Array of dimension sizes
            total = prod(obj.meshSizes);  % Total number of elements in the mesh
        
            % Create the grid indices for iteration
            [grid{1:numel(obj.breakpoints)}] = ndgrid(obj.breakpoints{:});
            
            % Flatten the grids into vectors for indexing in parfor
            gridVectors = cellfun(@(x) x(:), grid, 'UniformOutput', false);
        
            % Preallocate the mesh cell array
            mymesh = cell(obj.meshSizes);
            % Use parfor to populate the mesh
            fprintf('Calculating %.0f matrices \n', total);
            for i = 1:100:total
                fprintf("□");
            end
                fprintf("\n\n")

            parfor idx = 1:total
                % Extract the linearization values (linVals) for this index
                linVals = cellfun(@(vec) vec(idx), gridVectors, 'UniformOutput', false);
        
                % Compute the linearized system and store in the mesh
                sys = obj.getSS(linVals);
                sys = obj.getL(sys, 5);
                sys = obj.getX(sys, obj.Q, obj.R);
        
                % Assign to the mesh
                mymesh{idx} = sys;
        
                % Progress feedback (minimal overhead in parfor)
                if mod(idx, 100) == 0 || idx == total
                    fprintf('\b■\n', idx);
                end
            end
        
            % Assign to the object and save the result
            obj.mesh = mymesh;
            save("mesh.mat", "mymesh");
        end

    end

    methods (Access = private)
        
        function [linEq, linNum] = performLin(obj, eq, linVals)
            %PERFORMLINEARIZATION Helper to linearize equations
            %   eq: Array of symbolic equations to linearize
            %   subsVars: Cell array of variables to linearize
            %   subsVals: Struct with values for linearization

            % Substitute the linearization points

            % Linearize using Taylor expansion
            linEq = taylor(eq, obj.linVars, 'ExpansionPoint', linVals, 'Order', 2);

            % Substitute parameter values into the linearized model
            paramNames = fieldnames(obj.params);
            paramValues = struct2cell(obj.params);

            linNum = subs(linEq, paramNames, paramValues);
        end

    end
end



