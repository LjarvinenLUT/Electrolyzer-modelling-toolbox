clear;
close all;
clc;


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
weights = 'default';
switch i

    case 1 % Created test data
        
%         Uforfit = func.add(func.add(Uocv,Uact),func.add(Uohm,Ucon));
        
        [F,R,n_e] = getConstants();
        f = R/(F*n_e);
        alpha = 0.6; %
        j0 = 1e-6; % A/cm^2, exchange current density
        r = 0.1; % Ohm, total resistance
        j_lim = 1.5; % A/cm^2, limiting current density
        
        Workspace = struct('Variables',struct('T',T,'pCat',pH2,'pAn',pO2),'Coefficients',struct('alpha',alpha,'j0',j0,'r',r,'j_lim',j_lim));
        
        [SynData,FullData] = createSyntheticUI('workspace',Workspace,'jLims',[0.005 Inf],'jSigma',0.2,'uSigma',0.06);
        
%         Umeas = []; % "Measured" voltage
%         jmeas = []; % "Measured" current based on Buttler-Volmer equation
%         Tmeas = []; % "Measured" temperature with a drift
%         
%         T1 = T + (0:0.001:100)'*0;
%         U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...50, Activation overpotential sweep limits
%         Uocv = nernst(T1,pH2,pO2,'type',type);
%         
%         UohmFuncHandle = matlabFunction(str2sym(Uohm.equation));
%         UconFuncHandle = matlabFunction(str2sym(Ucon.equation));
%         U0 = Uocv.calculate;
%         
%         for ii = 1:length(U1)
%             jmeastemp = j0*(exp(alpha/(f*T1(ii))*U1(ii))-exp((alpha-1)/(f*T1(ii))*U1(ii))); % "Measured" current based on Buttler-Volmer equation
%             if jmeastemp >= j_lim
%                 break;
%             elseif jmeastemp > 5e-4
%                 U2 = UohmFuncHandle(jmeastemp,r); % Ohmic overpotential
%                 U3 = UconFuncHandle(F,R,T1(ii),jmeastemp,j_lim,n_e); % Concentration overpotential
%                 Umeastemp = U0(ii)+U1(ii)+U2+U3; % "Measured" voltage
%                 jmeas = [jmeas;jmeastemp];
%                 Umeas = [Umeas;Umeastemp];
%                 Tmeas = [Tmeas;T1(ii)];
%             end
%         end
%         pH2meas = pH2*ones(size(Tmeas));
%         pO2meas = pO2*ones(size(Tmeas));
%         
%         % Take samples from the dense data vectors
%         N = 20; % Number of evenly taken voltage samples
%         Usamples = linspace(min(Umeas)+0.2,max(Umeas),N)';
% %         jsamples = linspace(min(jmeas),jL-0.01,N)'; % Excluding mass transport limitations
% %         jmeassamp = nan(N,1);
% %         Umeassamp = nan(N,1);
%         iii = 1;
%         for ii = 1:N
%             Udif = abs(Umeas-Usamples(ii));
%             [~,ind] = min(Udif);
%             if iii == 1 || (iii > 1 && abs(jmeas(ind)-jmeassamp(iii-1))>0.02)
%                 jmeassamp(iii,1) = jmeas(ind); % Final sampled current vector
%                 Umeassamp(iii,1) = Umeas(ind); % Final sampled voltage vector
%                 Tmeassamp(iii,1) = Tmeas(ind); % Final sampled temperature vector with drift
%                 iii = iii+1;
%             end
%         end
%         
%         
%         % Adding sigma = 0.01*max normal error to measurements (Normal error seems to break the fit...)
%         p = 0.02;
%         jmeassamper = nan(length(jmeassamp),2);
%         Umeassamper = nan(length(jmeassamp),2);
%         for ii = 1:length(jmeassamp)
%             jm = jmeassamp(ii) + randn(200,1).*p*max(jmeassamp);
%             Um = Umeassamp(ii) + randn(200,1).*0.5*p*min(Umeassamp);
%             jmeassamper(ii,:) = [mean(jm) std(jm)];
%             Umeassamper(ii,:) = [mean(Um) std(Um)];
%         end
%         Uforfit.replaceParams('T',Tmeassamp);
        
        figure
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        hold on
        plot(FullData.current,FullData.voltage)
        title('Test "measurement" UI')
        xlabel('j')
        ylabel('U')
        
        %% Create electrolyzer model objects
        Emodel = electrolyzerModel('type',type); % Electrolyzer model containing the sampled temperature
        Emodel.setParams(struct('Variables',struct('T',T,'pCat',pH2,'pAn',pO2)));
%         Emodelfull = Emodel.copy; % Electrolyzer model containing the full dataset for temperature.
%         Emodelfull.replaceParams('T',T);
        Emodel.addPotentials('ocv','act','ohm','con');
%         Emodelfull.addPotentials('ocv','act','ohm','con');
        EmodelPS = Emodel.copy;
%         EmodelfullPS = Emodelfull.copy;
        
        %% Fit
        
        % Non-linear least squares error
        tic;
        [fitParam1,gof1] = Emodel.fitUI(SynData.voltage(:,1),SynData.current(:,1),'method','nllse','weights',weights);
        toc
		
		j0_fit1 = fitParam1.j0;
        alpha_fit1 = fitParam1.alpha;
        r_fit1 = fitParam1.r;
        j_lim_fit1 = fitParam1.j_lim;

		% Calculate voltage with calculate-method of func-object
%         Emodelfull.replaceParams(Emodel.potentialFunc.Workspace.Coefficients)
		Ufit1 = Emodel.calculate('current',FullData.current);
        
        RMSE1 = gof1.rmse; % Root mean squares error
        
        % Plotting
        
        figure
        hold on;
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        plot(FullData.current,Ufit1)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Non-linear least squares error')
        hold off;
        
        
        % Particleswarm
        tic;
        [fitParam2,gof2] = EmodelPS.fitUI(SynData.voltage(:,1),SynData.current(:,1),'method','ps','weights',weights);
        toc
        
        
		j0_fit2 = fitParam2.j0;
        alpha_fit2 = fitParam2.alpha;
        r_fit2 = fitParam2.r;
        j_lim_fit2 = fitParam2.j_lim;
        
		% Calculate voltage with calculate-method of func-object
%         EmodelfullPS.replaceParams(EmodelPS.potentialFunc.Workspace.Coefficients)
		Ufit2 = EmodelPS.calculate('current',FullData.current);
        
        RMSE2 = gof2.rmse; % Root mean squares error
        
        % Plotting
        
        figure
        hold on;
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        plot(FullData.current,Ufit2)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Particleswarm')
        hold off;
        
        %% Test with original parameters
%         Emodelfull2 = Emodelfull.copy;
%         Emodelfull2.replaceParams('alpha',alpha,'j0',j0,'r',r,'j_lim',j_lim)
        Emodel2 = Emodel.copy;
        Emodel2.replaceParams('alpha',alpha,'j0',j0,'r',r,'j_lim',j_lim)

        
        Ufit = Emodel2.calculate('current',FullData.current);
        
        RMSE = sqrt(mean((Emodel2.calculate('current',SynData.current(:,1))-SynData.voltage(:,1)).^2)); % Root mean squares error
        
        % Plotting
        
        figure
        hold on;
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        plot(FullData.current,Ufit)
        xlabel("j (A)")
        ylabel("U (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Optimal parameters')
        hold off;


        
    case 2 % Data from Jülich
        
        I0fit = nan(4,3);
        alphafit = nan(4,3);
        rfit = nan(4,3);
        RMSE = nan(2,3);
        Rsqrd = nan(2,3);
        
        load('JulichData.mat')
        
        %% Create electrolyzer model objects
        Emodel = electrolyzerModel('type','PEM');
        Emodel.setParams(struct('Variables',struct('T',T,'pCat',pH2,'pAn',pO2)));
        Emodel.addPotentials('ocv','act','ohm');
        
        for j = 1:3
            
            [Imeas,sortInd] = sort(JulichUI(j).I);
            I = [Imeas JulichUI(j).Istd(sortInd)];
            U = [JulichUI(j).U(sortInd) JulichUI(j).Ustd(sortInd)];
            
            %% Fit
            
            % Non-linear least squares error
            tic;
            [fitParam1,gof1] = Emodel.fitUI(U,I,'method','nllse','weights',weights);
            toc
            
            I0fit(1:2,j) = fitParam1.j0;
            alphafit(1:2,j) = fitParam1.alpha;
            rfit(1:2,j) = fitParam1.r;
            
            % Calculate voltage with calculate-method of func-object
            Ufit1 = Emodel.calculate('current',Imeas);
            
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
            [fitParam2,gof2] = Emodel.fitUI(U,I,'method','ps','weights',weights);
            toc
            
            
            I0fit(3:4,j) = fitParam2.j0;
            alphafit(3:4,j) = fitParam2.alpha;
            rfit(3:4,j) = fitParam2.r;
            
			% Calculate voltage with calculate-method of func-object
            Ufit2 = Emodel.calculate('current',Imeas);
            
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