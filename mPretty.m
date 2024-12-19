function mPretty(M)
    % Display matrix or array of matrices nicely formatted in the command window
    % Supports both symbolic and numeric matrices
    % If M is an array of matrices (e.g., [A, B, C, D]), they are displayed stacked vertically
    
    if isa(M, 'cell')
        % If M is a cell array of matrices, process each matrix in order
        for idx = 1:numel(M)
            fprintf('Matrix %d:\n', idx);
            mPretty(M{idx}); % Recursive call for each matrix
            fprintf('\n');    % Add spacing between matrices
        end
    else
        % Process a single matrix
        if isa(M, 'sym')
            % If M is symbolic, use the 'pretty' function
            pretty(M);
        else
            % If M is numeric, format and display each element
            [rows, cols] = size(M);
            for i = 1:rows
                fprintf('%8.4f ', M(i, :)); % Adjust width and precision as needed
                fprintf('\n');
            end
        end
    end
end
