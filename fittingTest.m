clear;
close all;
clc;

%% Global parameters
[F,R,n_e] = get_constants();

f = R/(n_e*F);

%%

for i = 4
    switch i
        case 1 % fit both activity and ohmic
            
            Uact = activation('model',1);
            
            Uohm = ohmic();
            
            Uforfit = @(j0,a,r,T,x) Uact(j0,a,T,x) + Uohm(r,x);
            
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[1e-10,0,0],...
                'Upper',[1,1,Inf],...
                'StartPoint',[0.001 0.5 1]);
            
            Ufitfun = fittype(Uforfit,...
                'dependent','z',...
                'coefficients',{'j0','a','r'},...
                'independent','x',...
                'problem','T',...
                'options',fo);
            
            %% Creating test data
            
            Tset = 289;
            U1 = ((3:0.01:10)*(f*Tset))'; % Uact/(f*T) = 0...10
            T = ones(size(U1))*Tset + rand(size(U1)); % T
            alpha = 0.4;
            j0 = 0.0001;
            x = j0*(exp(alpha./(f*T).*U1)-exp((alpha-1)./(f*T).*U1));
            r = 5;
            U2 = r*x;
            z = U1+U2;
            
            %% Fitting
            
            [fitted_curve,gof] = fit(x,z,Ufitfun,'problem',T);
            
            coeffvals = coeffvalues(fitted_curve);
            
            j0fit = coeffvals(1);
            afit = coeffvals(2);
            rfit = coeffvals(3);
            
            %% Plotting
            figure
            hold on;
            plot(x,z)
            plot(x,Uforfit(j0fit,afit,rfit,T,x))
            xlabel("j (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            hold off;
            
        case 2 % fit activity only
            
            Uact = activation('model',2);
            
            Uforfit = @(j0,a,T,x) Uact(j0,a,T,x);
            
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',[0,0],...
                'Upper',[1,1],...
                'StartPoint',[0.001 0.5]);
            
            Ufitfun = fittype(Uforfit,...
                'dependent','z',...
                'coefficients',{'j0','a'},...
                'independent','x',...
                'problem','T',...
                'options',fo);
            
            %% Creating test data
            
            Tset = 289;
            U1 = ((0:0.01:10)*(f*Tset))'; % Uact/(f*T) = 0...10
            T = ones(size(U1))*Tset + rand(size(U1)); % T
            alpha = 0.4;
            j0 = 0.0001;
            x = j0*(exp(alpha./(f*T).*U1)-exp((alpha-1)./(f*T).*U1));
            z = U1;
            
            %%
            
            [fitted_curve,gof] = fit(x,z,Ufitfun,'problem',T);
            
            coeffvals = coeffvalues(fitted_curve);
            
            j0fit = coeffvals(1);
            afit = coeffvals(2);
            
            %%
            figure
            hold on;
            plot(x,z)
            plot(x,Uforfit(j0fit,afit,T,x))
            xlabel("j (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            hold off;
            
            
            
            
        case 3 % fit ohmic only
            
            Uohm = @(r,j) r.*j;
            
            Uforfit = @(r,current) Uohm(r,current);
            
            Ufitfun = fittype(Uforfit,...
                'dependent','U',...
                'coefficients','r',...
                'independent','current'); % i and j are not accepted for some reason
            
            %% Creating test data
            
            T = 289;
            U1 = ((0:0.01:10)*(f*T))'; % Uact/(f*T) = 0...10
            alpha = 0.5;
            j0 = 0.0001;
            current = j0*(exp(alpha/(f*T).*U1)-exp((alpha-1)/(f*T).*U1));
            r = 5;
            U2 = current*r;
            U = U2;
            
            %% Fitting
            
            [fitted_curve,gof] = fit(current,U,Ufitfun);
            
            coeffvals = coeffvalues(fitted_curve);
            
            rfit = coeffvals(1);
            
            %% Plotting
            figure
            hold on;
            plot(current,U)
            plot(current,Uforfit(rfit,current))
            xlabel("j (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            hold off;
           
            
        case 4
            %% Creating function handle for fit
            
            T = 273.15+25; % Temperature
            type = "PEM"; % Cell type
            p1 = 50; % Hydrogen pressure
            p2 = 2; % Oxygen pressure
            Uocv = nernst(T,p1,p2,'type',type);
            Uact = activation(T);
            Uohm = ohmic();
            Ucon = concentration(T);
            
            Uforfit = @(j0,a,r,jL,Uerr,Current) Uocv + Uact(j0,a,Current) + Uohm(r,Current) + Ucon(jL,Current) + Uerr;
            
            %% Creating test data
            
            U1 = ((0.01:0.01:32)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential
            Tset = ones(size(U1))*T; % T vector
            alpha = 0.4;
            j0 = 1.3e-6;
            x = j0*(exp(alpha./(f*Tset).*U1)-exp((alpha-1)./(f*Tset).*U1)); % "Measured" current based on Buttler-Volmer equation
            r = 0.3;
            U2 = r*x; % Ohmic overpotential
            jL = 1.5;
            U3 = real(Ucon(jL,x));
            Uerr = 0.1;
            U4 = ones(size(U1))*Uerr; % Constant potential error
            U0 = ones(size(U1))*Uocv; % Open circuit voltage
            z = U0+U1+U2+U3+U4; % "Measured" voltage
            
            %% Fit
            
            fit_param = fit_UI(Uforfit,z,x);
            
            j0fit = fit_param(1);
            afit = fit_param(2);
            rfit = fit_param(3);
            jLfit = fit_param(4);
            Uerrfit = fit_param(5);
            
            Ufit = Uforfit(j0fit,afit,rfit,jLfit,Uerrfit,x);
            
            %% Plotting
            
            figure
            hold on;
            plot(x,z)
            plot(x,Ufit)
            xlabel("j (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            hold off;
    end
end

