% Concentration overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article

function Ucon = concentration(varargin)
    
    defaultModel = 1;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    
    %% Error checking
    
    
    %%
    
    switch model
        case 1
        case 2
        case 3
    end
end