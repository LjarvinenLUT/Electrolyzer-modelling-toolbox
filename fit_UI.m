

function fit_param = fit_UI(func_handle,Voltage,Current)

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[1e-10,0,0,0.01,-Inf],...
    'Upper',[1,1,Inf,Inf,Inf],...
    'StartPoint',[0.001 0.5 1 1 0]);

fitfun = fittype(func_handle,...
    'dependent','Voltage',...
    'coefficients',{'j0','a','r','jL','Uerr'},...
    'independent','Current',...
    'options',fo);

[fitted_curve,gof] = fit(Current,Voltage,fitfun);

fit_param = coeffvalues(fitted_curve);

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