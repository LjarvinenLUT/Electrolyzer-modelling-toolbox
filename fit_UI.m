function fit_param = fit_UI(func_handle,Voltage,Current)

coefficients = getFunctionArguments(func_handle);
[lower, upper, start] = getArgumentLimits(coefficients, Voltage, Current);

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',lower,...
    'Upper',upper,...
    'StartPoint',start);

fitfun = fittype(func_handle,...
    'dependent','Voltage',...
    'coefficients',coefficients,...
    'independent','Current',...
    'options',fo);

[fitted_curve,gof] = fit(Current,Voltage,fitfun);

coeff_values = coeffvalues(fitted_curve);

% Use map to store fitting parameters can be referenced by name
% Example: fit_param('j0')
fit_param = containers.Map(coefficients, coeff_values);
end


