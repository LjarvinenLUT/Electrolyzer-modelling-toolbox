

function fit_param = fit_UI(func_handle,Voltage,Current,varargin)

defaultMethod = 'PS';

parser = inputParser;
addRequired(parser,'func_handle',@(x) isa(x,'function_handle'))
addRequired(parser,'Voltage',@(x) isnumeric(x))
addRequired(parser,'Current',@(x) isnumeric(x))
addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x))

parse(parser,func_handle,Voltage,Current,varargin{:});

method = upper(string(parser.Results.method));

switch method
    
    case "NLLSE" % Non-Linear Least Squares Error regression approach
        % weight beginning and end
        x = (0:length(Current)-1)';
        weights = exp(x).^(-5/length(x)) + exp(x-x(end)).^(5/length(x));
        
        fo = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[1e-10,0,0,max(Current)],...
            'Upper',[1,1,Inf,Inf],...
            'StartPoint',[1e-6 0.5 1 max(Current)+1],...
            'Weights',weights,...
            'MaxIter',1500,...
            'Display','notify');
        
        fitfun = fittype(func_handle,...
            'dependent','Voltage',...
            'coefficients',{'j0','alpha','r','jL'},...
            'independent','Current',...
            'options',fo);
        
        [fitted_curve,gof] = fit(Current,Voltage,fitfun);
        
        fit_param = coeffvalues(fitted_curve);
        
      
    case "PS" % Particle swarm approach

        
        fun = @(x) sum((func_handle(x(1),x(2),x(3),x(4),Current)-Voltage).^2);
        nvars = 4;
        lb = [1e-10,0,0,max(Current)];
        ub = [1,1,40,5];
        
        options = optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon);
        
        [fit_param,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);
        
end

end


% numArguments = nargin( theFcn );
% 
% % Get the string description of the function
% functionString = func2str( theFcn );
% 
% % Allocate space for the cell-string
% arguments = cell( numArguments, 1 );
% 
% % The plan is to move a pair of indices along the "function string" field
% % looking for the commas. The names we want will be between these indices.
% %
% % We know from the form of the string that the first name starts two
% % characters after the "@"
% indexOfAtSign = find( functionString == '@', 1, 'first' );
% ai = indexOfAtSign + 2;
% % Therefore the first comma must be no sooner that the fourth character
% bi = 4;
% % When we start we have found no arguments
% numFound = 0;
% % We will keep looping until we have found all the arguments we expect
% while numFound < numArguments
%     % If we have found the end of a argument name
%     if functionString(bi) == ',' || functionString(bi) == ')'
%         % then increment the "numFound" counter
%         numFound = numFound+1;
%         % ... store the name
%         arguments{numFound} = functionString(ai:(bi-1));
%         % and increment the start index
%         ai = bi+1;
%         % Since the end must be beyond the start, we set the end index
%         % beyond the start index.
%         bi = ai+1;
%     else
%         % Otherwise increment the end index
%         bi = bi+1;
%     end
% end