clear;
close all;
clc;

addpath Utils

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



i = 3;
weights = 'none';
switch i

    case 1 
        %% Created test data
        
        Uforfit = combineFuncHandles({Uact,Uohm,Ucon,Uocv});
        
        alpha = 0.3; %
        j0 = 1e-6; % A/cm^2, exchange current density
        r = 0.1; % Ohm, total resistance
        j_lim = 1.5; % A/cm^2, limiting current density
        Uerr = 0; % V, constant voltage error
        
        Umeas = []; % "Measured" voltage
        jmeas = []; % "Measured" current based on Buttler-Volmer equation
        
        U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits
        
        for ii = 1:length(U1)
            jmeastemp = j0*(exp(alpha/(f*T)*U1(ii))-exp((alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
            if jmeastemp >= j_lim
                break;
            elseif jmeastemp > 1e-5
                U2 = Uohm(r,jmeastemp); % Ohmic overpotential
                U3 = Ucon(j_lim,jmeastemp); % Concentration overpotential
                Umeastemp = Uocv+U1(ii)+U2+U3+Uerr; % "Measured" voltage
                jmeas = [jmeas;jmeastemp];
                Umeas = [Umeas;Umeastemp];
            end
        end
        
        % Take samples from the dense data vectors
        N = 50; % Number of evenly taken current samples
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
        [fit_param1,fit_err1,gof1] = fit_UI(Uforfit,Umeassamper,jmeassamper,'method','nllse','weights',weights);
        toc
		
		j0fit1 = fit_param1.j0;
        alphafit1 = fit_param1.alpha;
        rfit1 = fit_param1.r;
        jLfit1 = fit_param1.j_lim;
        
		% Convert table to cell array which can then be used in a function call
		% instead of individual parameters
		fit_param1_cells = table2cell(fit_param1);

		% Use fit param cell array instead of individual parameters when calling
		% for voltage function
		Ufit1 = Uforfit(fit_param1_cells{:},jmeas);
        
        RMSE1 = gof1.rmse; % Root mean squares error
        
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
        [fit_param2,fit_err2,gof2] = fit_UI(Uforfit,Umeassamper,jmeassamper,'method','ps','weights',weights);
        toc
        
        
        j0fit2 = fit_param2.j0;
        alphafit2 = fit_param2.alpha;
        rfit2 = fit_param2.r;
        jLfit2 = fit_param2.j_lim;
        
		% Convert table to cell array which can then be used in a function call
		% instead of individual parameters
		fit_param2_cells = table2cell(fit_param2);

		% Use fit param cell array instead of individual parameters when calling
		% for voltage function
		Ufit2 = Uforfit(fit_param2_cells{:},jmeas);
        
        RMSE2 = gof2.rmse; % Root mean squares error
        
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
        
        Ufit = Uforfit(alpha,j0,j_lim,r,jmeas);
        
        RMSE = sqrt(mean((Uforfit(alpha,j0,j_lim,r,jmeassamp)-Umeassamper).^2)); % Root mean squares error
        
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


        
    case 2 
        %% Data from Jülich
        
        Uforfit = combineFuncHandles({Uocv,Uact,Uohm});
        
        I0fit = nan(2,3);
        alphafit = nan(2,3);
        rfit = nan(2,3);
        RMSE = nan(2,3);
        Rsqrd = nan(2,3);
        
        load('JulichData.mat')
        
        for j = 1:3
            
            [Imeas,sortInd] = sort(JulichUI(j).I);
            I = [Imeas JulichUI(j).Istd(sortInd)];
            U = [JulichUI(j).U(sortInd) JulichUI(j).Ustd(sortInd)];

            
            
            %% Fit
            
            % Non-linear least squares error
            tic;
            [fit_param1,fit_err1,gof1] = fit_UI(Uforfit,U,I,'method','nllse','weights',weights);
            toc
            
            I0fit(1,j) = fit_param1.j0;
            alphafit(1,j) = fit_param1.alpha;
            rfit(1,j) = fit_param1.r;
            
			% Convert table to cell array which can then be used in a function call
			% instead of individual parameters
			fit_param1_cells = table2cell(fit_param1);

			% Use fit param cell array instead of individual parameters when calling
			% for voltage function
            Ufit1 = Uforfit(fit_param1_cells{:},I(:,1));
            
            RMSE(1,j) = gof1.rmse; % Root mean squares error
            Rsqrd(1,j) = gof1.rsquare; % R^2 of the fit
            
            % Plotting
            
            figure
            hold on;
            errorbar(I(:,1),U(:,1),U(:,2),U(:,2),I(:,2),I(:,2),'o')
            plot(I(:,1),Ufit1)
            xlabel("I (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Non-linear least squares error: Jülich ' num2str(j)])
            hold off;
            
            
            % Particleswarm
            tic;
            [fit_param2,fit_err2,gof2] = fit_UI(Uforfit,U(:,1),I(:,1),'method','ps','weights',weights);
            toc
            
            
            I0fit(2,j) = fit_param2.j0;
            alphafit(2,j) = fit_param2.alpha;
            rfit(2,j) = fit_param2.r;
            
			% Convert table to cell array which can then be used in a function call
			% instead of individual parameters
			fit_param2_cells = table2cell(fit_param2);

			% Use fit param cell array instead of individual parameters when calling
			% for voltage function
            Ufit2 = Uforfit(fit_param2_cells{:},I(:,1));
            
            RMSE(2,j) = gof2.rmse; % Root mean squares error
            Rsqrd(2,j) = gof2.rsquare; % R^2 of the fit
            
            % Plotting
            
            figure
            hold on;
            errorbar(I(:,1),U(:,1),U(:,2),U(:,2),I(:,2),I(:,2),'o')
            plot(Imeas,Ufit2)
            xlabel("I (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Particleswarm: Jülich ' num2str(j)])
            hold off;
            
        end
        
    
    case 3 
        %% Data from Baby Piel
        
        load('AlkaliData.mat')
        
        T = mean(AlkaliUI.T(:,1));
        P = mean(AlkaliUI.P(:,1));
        
        type = "alkaline"; % Cell type
        m = 7.64; % electrolyte molality for 30% KOH solution
        Uocv = nernst(T,P,m,'type',type,'electrolyte','NaOH');
        activation_model = 2;
        Uact = activation(T,'model',activation_model);
        Uohm = ohmic();
        
        Uforfit = combineFuncHandles({Uocv,Uact,Uohm});
        
        %% Fit
        
        % Non-linear least squares error
        tic;
        [fit_param1,fit_err1,gof1] = fit_UI(Uforfit,AlkaliUI.U,AlkaliUI.j,'method','nllse','weights',weights);
        toc
        
        I0fit(1,1) = fit_param1.j0;
        if activation_model == 1
            alphafit(1,1) = 1/2;
        else
            alphafit(1,1) = fit_param1.alpha;
        end
        rfit(1,1) = fit_param1.r;
        
        % Convert table to cell array which can then be used in a function call
        % instead of individual parameters
        fit_param1_cells = table2cell(fit_param1);
        
        % Use fit param cell array instead of individual parameters when calling
        % for voltage function
        test_j = 0.001:0.001:max(AlkaliUI.j(:,1));
        Ufit1 = Uforfit(fit_param1_cells{:},test_j);
        
        RMSE(1,1) = gof1.rmse; % Root mean squares error
        Rsqrd(1,1) = gof1.rsquare; % R^2 of the fit
        
        % Plotting
        
        figure
        hold on;
        errorbar(AlkaliUI.j(:,1),AlkaliUI.U(:,1),AlkaliUI.U(:,2),AlkaliUI.U(:,2),AlkaliUI.j(:,2),AlkaliUI.j(:,2),'o')
        plot(test_j,Ufit1)
        xlabel("I (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Non-linear least squares error: Alkali')
        hold off;
        
        
        % Particleswarm
        tic;
        [fit_param2,fit_err2,gof2] = fit_UI(Uforfit,AlkaliUI.U(:,1),AlkaliUI.j(:,1),'method','ps','weights',weights);
        toc
        
        
        I0fit(2,1) = fit_param2.j0;
        if activation_model == 1
            alphafit(2,1) = 1/2;
        else
            alphafit(2,1) = fit_param2.alpha;
        end
        rfit(2,1) = fit_param2.r;
        
        % Convert table to cell array which can then be used in a function call
        % instead of individual parameters
        fit_param2_cells = table2cell(fit_param2);
        
        % Use fit param cell array instead of individual parameters when calling
        % for voltage function
        Ufit2 = Uforfit(fit_param2_cells{:},test_j);
        
        RMSE(2,1) = gof2.rmse; % Root mean squares error
        Rsqrd(2,1) = gof2.rsquare; % R^2 of the fit
        
        % Plotting
        
        figure
        hold on;
        errorbar(AlkaliUI.j(:,1),AlkaliUI.U(:,1),AlkaliUI.U(:,2),AlkaliUI.U(:,2),AlkaliUI.j(:,2),AlkaliUI.j(:,2),'o')
        plot(test_j,Ufit2)
        xlabel("I (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Particleswarm: Alkali')
        hold off;
        
end