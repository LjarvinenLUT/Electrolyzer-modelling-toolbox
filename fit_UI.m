

function fit_param = fit_UI(func_handle,Voltage,Current)

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[1e-10,0,0,-Inf],...
    'Upper',[1,1,Inf,Inf],...
    'StartPoint',[0.001 0.5 1 0]);

fitfun = fittype(func_handle,...
    'dependent','Voltage',...
    'coefficients',{'j0','a','r','Uerr'},...
    'independent','Current',...
    'options',fo);

[fitted_curve,gof] = fit(Current,Voltage,fitfun);

fit_param = coeffvalues(fitted_curve);

end