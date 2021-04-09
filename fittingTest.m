clear;
close all;
clc;

%% Global parameters
[F,R,n_e] = get_constants();

f = R/(n_e*F);


%% Creating function handle for fit

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pH2 = 30; % bar
pO2 = 2; % bar
Uocv = nernst(T,pH2,pO2,'type',type);
Uact = activation(T,'model',3);
Uohm = ohmic();
Ucon = concentration(T);

Uforfit = @(j0,alpha,r,jL,Current) Uocv + Uact(j0,alpha,Current) + Uohm(r,Current) + Ucon(jL,Current);

i = 2;
switch i

    case 1 % Created test data
        
        
        %             U1 = ((0.01:0.01:32)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential
        %             Tset = ones(size(U1))*T; % T vector
        %             alpha = 0.4;
        %             j0 = 1.3e-6;
        %             x = j0*(exp(alpha./(f*Tset).*U1)-exp((alpha-1)./(f*Tset).*U1)); % "Measured" current based on Buttler-Volmer equation
        %             r = 0.3;
        %             U2 = r*x; % Ohmic overpotential
        %             jL = 1.5;
        %             U3 = real(Ucon(jL,x));
        %
        %             U4 = ones(size(U1))*Uerr; % Constant potential error
        %             U0 = ones(size(U1))*Uocv; % Open circuit voltage
        %             z = U0+U1+U2+U3+U4; % "Measured" voltage
        
        alpha = 0.5; %
        j0 = 1e-6; % A/cm^2, exchange current density
        r = 0.3; % Ohm, total resistance
        jL = 1.5; % A/cm^2, limiting current density
        Uerr = 0; % V, constant voltage error
        
        Umeas = []; % "Measured" voltage
        jmeas = []; % "Measured" current based on Buttler-Volmer equation
        
        U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits
        
        for ii = 1:length(U1)
            jmeastemp = j0*(exp(alpha/(f*T)*U1(ii))-exp((alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
            if jmeastemp >= jL
                break;
            elseif jmeastemp > 0.001
                U2 = Uohm(r,jmeastemp); % Ohmic overpotential
                U3 = Ucon(jL,jmeastemp); % Concentration overpotential
                Umeastemp = Uocv+U1(ii)+U2+U3+Uerr; % "Measured" voltage
                jmeas = [jmeas;jmeastemp];
                Umeas = [Umeas;Umeastemp];
            end
        end
        
        % Take samples from the dense data vectors
        N = 100; % Number of evenly taken current samples
        % jsamples = linspace(min(jmeas),max(jmeas),N)';
        jsamples = linspace(min(jmeas),jL-0.01,N)'; % Excluding mass transport limitations
        jmeassamp = nan(N,1);
        Umeassamp = nan(N,1);
        for ii = 1:N
            jdif = abs(jmeas-jsamples(ii));
            [~,ind] = min(jdif);
            jmeassamp(ii) = jmeas(ind); % Final sampled current vector
            Umeassamp(ii) = Umeas(ind); % Final sampled voltage vector
        end
        
        
        % Adding p*100% error to measurements
        p = 0.01;
        jmeassamp = jmeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        Umeassamp = Umeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        
        
        figure
        scatter(jmeassamp,Umeassamp)
        hold on
        plot(jmeas,Umeas)
        title('Test "measurement" UI')
        xlabel('j')
        ylabel('U')
        
        %% Fit
        
        % Non-linear least squares error
        tic;
        fit_param1 = fit_UI(Uforfit,Umeassamp,jmeassamp,'method','nllse');
        toc
        
        
        j0fit1 = fit_param1('j0');
        alphafit1 = fit_param1('alpha');
        rfit1 = fit_param1('r');
        jLfit1 = fit_param1('jL');
        % j0fit = j0;
        % alphafit = alpha;
        % rfit = r;
        % jLfit = jL;
        
        Ufit1 = Uforfit(j0fit1,alphafit1,rfit1,jLfit1,jmeas);
        
        MSE1 = mean((Ufit1-Umeas).^2); % Mean squares error
        
        % Plotting
        
        figure
        hold on;
        scatter(jmeassamp,Umeassamp)
        plot(jmeas,Ufit1)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Non-linear least squares error')
        hold off;
        
        
        % Particleswarm
        tic;
        fit_param2 = fit_UI(Uforfit,Umeassamp,jmeassamp,'method','ps');
        toc
        
        
        j0fit2 = fit_param2('j0');
        alphafit2 = fit_param2('alpha');
        rfit2 = fit_param2('r');
        jLfit2 = fit_param2('jL');
        % j0fit = j0;
        % alphafit = alpha;
        % rfit = r;
        % jLfit = jL;
        
        Ufit2 = Uforfit(j0fit2,alphafit2,rfit2,jLfit2,jmeas);
        
        MSE2 = mean((Ufit2-Umeas).^2); % Mean squares error
        
        % Plotting
        
        figure
        hold on;
        scatter(jmeassamp,Umeassamp)
        plot(jmeas,Ufit2)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Particleswarm')
        hold off;
        
        %% Test with original parameters
        
        j0fit = j0;
        alphafit = alpha;
        rfit = r;
        jLfit = jL;
        
        Ufit = Uforfit(j0fit,alphafit,rfit,jLfit,jmeas);
        
        MSE = mean((Ufit-Umeas).^2); % Mean squares error
        
        % Plotting
        
        figure
        hold on;
        scatter(jmeassamp,Umeassamp)
        plot(jmeas,Ufit)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Optimal parameters')
        hold off;


        
    case 2 % Data from Jülich
        
        j0fit = nan(2,3);
        alphafit = nan(2,3);
        jLfit = nan(2,3);
        rfit = nan(2,3);
        MSE = nan(2,3);
        
        load('JulichData.mat')
        
        for j = 1:3
            
            [jmeas,sortInd] = sort(JulichUI(j).I);
            Umeas = JulichUI(j).U(sortInd);
            
            
            %% Fit
            
            % Non-linear least squares error
            tic;
            fit_param1 = fit_UI(Uforfit,Umeas,jmeas,'method','nllse');
            toc
            
            
            j0fit1 = fit_param1('j0');
            alphafit1 = fit_param1('alpha');
            rfit1 = fit_param1('r');
            jLfit1 = fit_param1('jL');
            
            j0fit(1,j) = j0fit1;
            alphafit(1,j) = alphafit1;
            rfit(1,j) = rfit1;
            jLfit(1,j) = jLfit1;
            
            Ufit1 = Uforfit(j0fit1,alphafit1,rfit1,jLfit1,jmeas);
            
            MSE1 = mean((Ufit1-Umeas).^2); % Mean squares error
            MSE(1,j) = MSE1;
            
            % Plotting
            
            figure
            hold on;
            scatter(jmeas,Umeas)
            plot(jmeas,Ufit1)
            xlabel("j (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Non-linear least squares error: Jülich ' num2str(j)])
            hold off;
            
            
            % Particleswarm
            tic;
            fit_param2 = fit_UI(Uforfit,Umeas,jmeas,'method','ps');
            toc
            
            
            j0fit2 = fit_param2('j0');
            alphafit2 = fit_param2('alpha');
            rfit2 = fit_param2('r');
            jLfit2 = fit_param2('jL');
            
            j0fit(2,j) = j0fit2;
            alphafit(2,j) = alphafit2;
            rfit(2,j) = rfit2;
            jLfit(2,j) = jLfit2;
            
            Ufit2 = Uforfit(j0fit2,alphafit2,rfit2,jLfit2,jmeas);
            
            MSE2 = mean((Ufit2-Umeas).^2); % Mean squares error
            MSE(2,j) = MSE2;
            
            % Plotting
            
            figure
            hold on;
            scatter(jmeas,Umeas)
            plot(jmeas,Ufit2)
            xlabel("j (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Particleswarm: Jülich ' num2str(j)])
            hold off;
            
        end
end

