clear;
close all;
clc;

%% Global parameters
N_A = 6.02214076e23; % 1/mol, Avogadro constant
k_B = 1.380649e-23; % J/K, Boltzmann constant
Q_e = 1.602176634e-19; % C, Elemental charge
global F n_e R
F = Q_e*N_A; % C/mol, Faraday constant
R = k_B*N_A; % J/K/mol, Universal gas constant
n_e = 2; % Number of electrons transferred in one reaction

f = R/(n_e*F);

%%

i = 3;

switch i
    case 1 % fit both activity and ohmic
        
        Uact = activation('model',1);
        
        Uohm = @(r,j) r.*j;
        
        Uforfit = @(j0,a,r,T,x) Uact(j0,a,T,x) + Uohm(r,x);
        
        fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
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
        U1 = ((0:0.01:10)*(f*Tset))'; % Uact/(f*T) = 0...10
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
        
        Uact = activation();
        
        y = 289;
        
        Uforfit = @(j0,T,x) Uact(j0,T,x);
        
        Ufitfun = fittype(Uforfit,...
            'dependent','y',...
            'coefficients','j0',...
            'independent','x',...
            'problem','T');
        
        %% Creating test data
        
        U1 = ((0:0.01:10)*(f*y))'; % Uact/(f*T) = 0...10
        alpha = 0.5;
        j0 = 0.0001;
        x = j0*(exp(alpha/(f*y).*U1)-exp((alpha-1)/(f*y).*U1));
        y = U1;
        
        %%
        
        [fitted_curve,gof] = fit(x,y,Ufitfun,'problem',y);
        
        coeffvals = coeffvalues(fitted_curve);
        
        j0fit = coeffvals(1);
        
        %%
        figure
        hold on;
        plot(x,y)
        plot(x,Uforfit(j0fit,y,x))
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
        
        y = 289;
        U1 = ((0:0.01:10)*(f*y))'; % Uact/(f*T) = 0...10
        alpha = 0.5;
        j0 = 0.0001;
        current = j0*(exp(alpha/(f*y).*U1)-exp((alpha-1)/(f*y).*U1));
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
        
end

