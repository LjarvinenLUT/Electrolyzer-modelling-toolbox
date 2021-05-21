clear;
close all;
clc;

%% Global parameters
[F,R,n_e] = getConstants();

f = R/(n_e*F);


%% Creating function handle for fit

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pH2 = 30; % bar
pO2 = 2; % bar
Uocv = nernst(T,pH2,pO2,'type',type);
Uact = activation('model',2);
Uohm = ohmic();
Ucon = concentration();



i = 1;
weights = 'none';
switch i

    case 1 % Created test data
        
        Uforfit = addFuncs(addFuncs(Uocv,Uact),addFuncs(Uohm,Ucon));
        
        
        alpha = 0.3; %
        j0 = 1e-6; % A/cm^2, exchange current density
        r = 0.1; % Ohm, total resistance
        j_lim = 1.5; % A/cm^2, limiting current density
        
        
        Umeas = []; % "Measured" voltage
        jmeas = []; % "Measured" current based on Buttler-Volmer equation
        
        U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits
        
        UohmFuncHandle = matlabFunction(str2sym(Uohm.equation));
        UconFuncHandle = matlabFunction(str2sym(Ucon.equation));
        U0 = Uocv.calculate;
        
        for ii = 1:length(U1)
            jmeastemp = j0*(exp(alpha/(f*T)*U1(ii))-exp((alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
            if jmeastemp >= j_lim
                break;
            elseif jmeastemp > 1e-5
                U2 = UohmFuncHandle(jmeastemp,r); % Ohmic overpotential
                U3 = UconFuncHandle(F,R,T,jmeastemp,j_lim,n_e); % Concentration overpotential
                Umeastemp = U0+U1(ii)+U2+U3; % "Measured" voltage
                jmeas = [jmeas;jmeastemp];
                Umeas = [Umeas;Umeastemp];
            end
        end
        Tmeas = T*ones(size(jmeas));
        pH2meas = Uforfit.Workspace.Variables.pH2.*ones(size(jmeas));
        pO2meas = Uforfit.Workspace.Variables.pO2.*ones(size(jmeas));
        
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
        Tmeassamp = T*ones(size(jmeassamp));
        pH2meassamp = Uforfit.Workspace.Variables.pH2.*ones(size(jmeassamp));
        pO2meassamp = Uforfit.Workspace.Variables.pO2.*ones(size(jmeassamp));
        
        
        % Adding p*100% error to measurements
        p = 0.01;
        jmeassamper = jmeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        Umeassamper = Umeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        Tmeassamper = Tmeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
        Uforfit.Workspace.Variables.T = Tmeassamper;
        Uforfit.Workspace.Variables.pH2 = pH2meassamp;
        Uforfit.Workspace.Variables.pO2 = pO2meassamp;
        
        figure
        scatter(jmeassamp,Umeassamp)
        hold on
        plot(jmeas,Umeas)
        title('Test "measurement" UI')
        xlabel('j')
        ylabel('U')
        
        %% Create electrolyzer model objects
        Emodel1 = electrolyzerModel('type',type);
        Emodel1.addPotential(Uforfit.copy);
        Emodel2 = Emodel1.copy;
        
        %% Fit
        
        % Non-linear least squares error
        tic;
        [fitParam1,gof1] = Emodel1.fitUI(Umeassamper,jmeassamper,'method','nllse','weights',weights);
        toc
		
		j0_fit1 = fitParam1.j0(1);
        alpha_fit1 = fitParam1.alpha(1);
        r_fit1 = fitParam1.r(1);
        j_lim_fit1 = fitParam1.j_lim(1);

		% Calculate voltage with calculate-method of func-object
		Ufit1 = Emodel1.calculate('current',jmeas,'T',Tmeas,'pH2',pH2meas,'pO2',pO2meas);
        
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
        [fitParam2,gof2] = Emodel2.fitUI(Umeassamper,jmeassamper,'method','ps','weights',weights);
        toc
        
        
		j0_fit2 = fitParam2.j0(1);
        alpha_fit2 = fitParam2.alpha(1);
        r_fit2 = fitParam2.r(1);
        j_lim_fit2 = fitParam2.j_lim(1);
        
		% Calculate voltage with calculate-method of func-object
		Ufit2 = Emodel2.calculate('current',jmeas,'T',Tmeas,'pH2',pH2meas,'pO2',pO2meas);
        
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
        
        Uforfit.Workspace.Coefficients.alpha = alpha ; %
        Uforfit.Workspace.Coefficients.j0 = j0; % A/cm^2, exchange current density
        Uforfit.Workspace.Coefficients.r = r; % Ohm, total resistance
        Uforfit.Workspace.Coefficients.j_lim = j_lim; % A/cm^2, limiting current density
        
        Ufit = Uforfit.calculate('current',jmeas,'T',Tmeas,'pH2',pH2meas,'pO2',pO2meas);
        
        RMSE = sqrt(mean((Uforfit.calculate('current',jmeassamp)-Umeassamper).^2)); % Root mean squares error
        
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
        
        Uforfit = addFuncs(addFuncs(Uocv,Uact),Uohm);
        
        
        
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

            Uforfit1 = Uforfit.copy;
            Uforfit2 = Uforfit.copy;
            
            %% Fit
            
            % Non-linear least squares error
            tic;
            [fitParam1,gof1] = fitUI(Uforfit1,U,I,'method','nllse','weights',weights);
            toc
            
            I0fit(1,j) = fitParam1.j0(1);
            alphafit(1,j) = fitParam1.alpha(1);
            rfit(1,j) = fitParam1.r(1);
            
            % Calculate voltage with calculate-method of func-object
            Ufit1 = Uforfit1.calculate('current',Imeas);
            
            RMSE(1,j) = gof1.rmse; % Root mean squares error
            Rsqrd(1,j) = gof1.rsquare; % R^2 of the fit
            
            % Plotting
            
            figure
            hold on;
            errorbar(I(:,1),U(:,1),U(:,2),U(:,2),I(:,2),I(:,2),'o')
            plot(Imeas,Ufit1)
            xlabel("I (A)")
            ylabel("U (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Non-linear least squares error: Jülich ' num2str(j)])
            hold off;
            
            
            % Particleswarm
            tic;
            [fitParam2,gof2] = fitUI(Uforfit2,U(:,1),I(:,1),'method','ps','weights',weights);
            toc
            
            
            I0fit(2,j) = fitParam2.j0(1);
            alphafit(2,j) = fitParam2.alpha(1);
            rfit(2,j) = fitParam2.r(1);
            
			% Calculate voltage with calculate-method of func-object
            Ufit2 = Uforfit2.calculate('current',Imeas);
            
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
end