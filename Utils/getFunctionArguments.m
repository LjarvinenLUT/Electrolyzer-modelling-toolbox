function [coefficients] = getFunctionArguments(theFcn)
    %GETFUNCTIONARGUMENTS returns the fucntion arguments from the provided
    %function as a cell array.
    
    % Get the amount of arguments
    numArguments = nargin( theFcn );

    % Get the string description of the function
    functionString = func2str( theFcn );

    % Allocate space for the cell-string. One of the arguments is the
    % independent fit variable which is not needed here
    coefficients = cell( 1, numArguments - 1);

    % The plan is to move a pair of indices along the "function string" field
    % looking for the commas. The names we want will be between these indices.
    %
    % We know from the form of the string that the first name starts two
    % characters after the "@"
    indexOfAtSign = find( functionString == '@', 1, 'first' );
    ai = indexOfAtSign + 2;
    % Therefore the first comma must be no sooner that the fourth character
    bi = 4;
    % When we start we have found no arguments
    numFound = 0;
    % We will keep looping until we have found all the arguments we expect
    while numFound < numArguments
        % If we have found the end of a argument name
        if functionString(bi) == ',' || functionString(bi) == ')'
            % then increment the "numFound" counter
            numFound = numFound+1;
            % ... store the name
            if ~strcmpi(functionString(ai:(bi-1)), "Current")
                coefficients{numFound} = functionString(ai:(bi-1));
            end
            % and increment the start index
            ai = bi+1;
            % Since the end must be beyond the start, we set the end index
            % beyond the start index.
            bi = ai+1;
        else
            % Otherwise increment the end index
            bi = bi+1;
        end
    end
end