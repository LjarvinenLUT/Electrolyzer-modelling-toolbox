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
Uact = activation(T,'model',2);
Uohm = ohmic();
Ucon = concentration(T);



i = 2;
switch i

    case 1 % Created test data
        
        Uforfit = @(j0,alpha,r,jL,Current) Uocv + Uact(j0,alpha,Current) + Uohm(r,Current) + Ucon(jL,Current);
        
        alpha = 0.5; %
        j0 = 1e-6; % A/cm^2, exchange current density
        r = 0.1; % Ohm, total resistance
        jL = 1.5; % A/cm^2, limiting current density
        Uerr = 0; % V, constant voltage error
        
        Umeas = []; % "Measured" voltage
        jmeas = []; % "Measured" current based on Buttler-Volmer equation
        
        U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits
        
        for ii = 1:length(U1)
            jmeastemp = j0*(exp(alpha/(f*T)*U1(ii))-exp((alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
            if jmeastemp >= jL
                break;
            elseif jmeastemp > 1e-5
                U2 = Uohm(r,jmeastemp); % Ohmic overpotential
                U3 = Ucon(jL,jmeastemp); % Concentration overpotential
                Umeastemp = Uocv+U1(ii)+U2+U3+Uerr; % "Measured" voltage
                jmeas = [jmeas;jmeastemp];
                Umeas = [Umeas;Umeastemp];
            end
        end
        
        % Take samples from the dense data vectors
        N = 25; % Number of evenly taken current samples
        jsamples = linspace(min(jmeas),max(jmeas),N)';
%         jsamples = linspace(min(jmeas),jL-0.01,N)'; % Excluding mass transport limitations
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
        jmeassamper = jmeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        Umeassamper = Umeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        
        
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
        [fit_param1,fit_err1] = fit_UI(Uforfit,Umeassamper,jmeassamper,'method','nllse');
        toc
		
		j0fit1 = fit_param1.j0;
        alphafit1 = fit_param1.alpha;
        rfit1 = fit_param1.r;
        jLfit1 = fit_param1.jL;
        
		% Convert table to cell array which can then be used in a function call
		% instead of individual parameters
		fit_param1_cells = table2cell(fit_param1);

		% Use fit param cell array instead of individual parameters when calling
		% for voltage function
		Ufit1 = Uforfit(fit_param1_cells{:},jmeas);
        
        RMSE1 = mean((Uforfit(fit_param1_cells{:},jmeassamp)-Umeassamper).^2); % Mean squares error
        
        % Plotting
        
        figure
        hold on;
        scatter(jmeassamper,Umeassamper)
        plot(jmeas,Ufit1)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Non-linear least squares error')
        hold off;
        
        
        % Particleswarm
        tic;
        [fit_param2,fit_err2] = fit_UI(Uforfit,Umeassamper,jmeassamper,'method','ps');
        toc
        
        
        j0fit2 = fit_param2.j0;
        alphafit2 = fit_param2.alpha;
        rfit2 = fit_param2.r;
        jLfit2 = fit_param2.jL;
        
		% Convert table to cell array which can then be used in a function call
		% instead of individual parameters
		fit_param2_cells = table2cell(fit_param2);

		% Use fit param cell array instead of individual parameters when calling
		% for voltage function
		Ufit2 = Uforfit(fit_param2_cells{:},jmeas);
        
        RMSE2 = mean((Uforfit(fit_param2_cells{:},jmeassamp)-Umeassamper).^2); % Mean squares error
        
        % Plotting
        
        figure
        hold on;
        scatter(jmeassamper,Umeassamper)
        plot(jmeas,Ufit2)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Particleswarm')
        hold off;
        
        %% Test with original parameters
        
        Ufit = Uforfit(j0,alpha,r,jL,jmeas);
        
        RMSE = mean((Uforfit(j0,alpha,r,jL,jmeassamp)-Umeassamper).^2); % Mean squares error
        
        % Plotting
        
        figure
        hold on;
        scatter(jmeassamper,Umeassamper)
        plot(jmeas,Ufit)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Optimal parameters')
        hold off;


        
    case 2 % Data from Jülich
        
        Uforfit = @(j0,alpha,r,Current) Uocv + Uact(j0,alpha,Current) + Uohm(r,Current);
        
        I0fit = nan(2,3);
        alphafit = nan(2,3);
        rfit = nan(2,3);
        RMSE = nan(2,3);
        
        load('JulichData.mat')
        
        for j = 1:3
            
            [Imeas,sortInd] = sort(JulichUI(j).I);
            Umeas = JulichUI(j).U(sortInd);
            
            
            %% Fit
            
            % Non-linear least squares error
            tic;
            [fit_param1,fit_err1,gof1] = fit_UI(Uforfit,Umeas,Imeas,'method','nllse');
            toc
            
            I0fit(1,j) = fit_param1.j0;
            alphafit(1,j) = fit_param1.alpha;
            rfit(1,j) = fit_param1.r;
            
			% Convert table to cell array which can then be used in a function call
			% instead of individual parameters
			fit_param1_cells = table2cell(fit_param1);

			% Use fit param cell array instead of individual parameters when calling
			% for voltage function
            Ufit1 = Uforfit(fit_param1_cells{:},Imeas);
            
            RMSE1 = gof1.rmse; % Root mean squares error
            RMSE(1,j) = RMSE1;
            
            % Plotting
            
            figure
            hold on;
            scatter(Imeas,Umeas)
            plot(Imeas,Ufit1)
            xlabel("I (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Non-linear least squares error: Jülich ' num2str(j)])
            hold off;
            
            
            % Particleswarm
            tic;
            [fit_param2,fit_err2,gof2] = fit_UI(Uforfit,Umeas,Imeas,'method','ps');
            toc
            
            
            I0fit(2,j) = fit_param2.j0;
            alphafit(2,j) = fit_param2.alpha;
            rfit(2,j) = fit_param2.r;
            
			% Convert table to cell array which can then be used in a function call
			% instead of individual parameters
			fit_param2_cells = table2cell(fit_param2);

			% Use fit param cell array instead of individual parameters when calling
			% for voltage function
            Ufit2 = Uforfit(fit_param2_cells{:},Imeas);
            
            RMSE2 = mean((Ufit2-Umeas).^2); % Mean squares error
            RMSE(2,j) = RMSE2;
            
            % Plotting
            
            figure
            hold on;
            scatter(Imeas,Umeas)
            plot(Imeas,Ufit2)
            xlabel("I (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Particleswarm: Jülich ' num2str(j)])
            hold off;
            
        end
end