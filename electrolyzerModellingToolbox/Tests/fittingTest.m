clear;
close all;
clc;


%% Creating function handle for fit

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pH2 = 30; % bar
pO2 = 2; % bar
Uocv = nernst(type);
Uact = activation('model',2);
Uohm = ohmic();
Ucon = concentration();

eModel = electrolyzerModel('type',type); % Electrolyzer model containing the sampled temperature
eModel.setInWorkspace(struct('Variables',struct('T',T,'pCat',pH2,'pAn',pO2)));
%%
i = 1;
switch i

    case 1 % Created test data
        
%         Uforfit = func.add(func.add(Uocv,Uact),func.add(Uohm,Ucon));
        
        [F,R,n_e] = getConstants();
        f = R/(F*n_e);
        alpha = 0.6; %
        j0 = 1e-6; % A/cm^2, exchange current density
        r = 0.1; % Ohm, total resistance
        j_lim = 1.5; % A/cm^2, limiting current density
        
        synModel = eModel.copy;
        synModel.setInWorkspace(struct('Parameters',struct('alpha',alpha,'j0',j0,'r',r,'j_lim',j_lim)))
        
        [SynData,FullData] = createSyntheticUI('model',synModel,'jLims',[0.01 j_lim-0.01],'jErr',0.01,'uErr',0.01);
        
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
%         Uforfit.replaceInWorkspace('T',Tmeassamp);
        
        figure
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        hold on
        plot(FullData.current,FullData.voltage)
        title('Test "measurement" UI')
        xlabel('j')
        ylabel('U')
        
        %% Create electrolyzer model objects
%         eModel = electrolyzerModel('type',type); % Electrolyzer model containing the sampled temperature
%         eModel.setInWorkspace(struct('Variables',struct('T',T,'pCat',pH2,'pAn',pO2)));
%         Emodelfull = Emodel.copy; % Electrolyzer model containing the full dataset for temperature.
%         Emodelfull.replaceInWorkspace('T',T);
        eModel.addFuncs('ocv','act','ohm','con');
%         Emodelfull.addFuncs('ocv','act','ohm','con');
        eModel2 = eModel.copy;
%         EmodelfullPS = Emodelfull.copy;
        
        %% Fit
        
        weights = 'hl';
        
        % Non-linear least squares error
        tic;
        [fitParams1,gof1] = eModel.fitUI(SynData.voltage(:,1),SynData.current(:,1),'method','nllse','weights',weights);
        toc
		
		j0_fit1 = fitParams1.j0;
        alpha_fit1 = fitParams1.alpha;
        r_fit1 = fitParams1.r;
        j_lim_fit1 = fitParams1.j_lim;

		% Calculate voltage with calculate-method of func-object
%         Emodelfull.replaceInWorkspace(Emodel.potentialFunc.Workspace.Parameters)
		Ufit1 = eModel.calculate('current',FullData.current);
        
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
        [fitParams2,gof2] = eModel2.fitUI(SynData.voltage(:,1),SynData.current(:,1),'method','ps','weights',weights);
        toc
        
        
		j0_fit2 = fitParams2.j0;
        alpha_fit2 = fitParams2.alpha;
        r_fit2 = fitParams2.r;
        j_lim_fit2 = fitParams2.j_lim;
        
		% Calculate voltage with calculate-method of func-object
%         EmodelfullPS.replaceInWorkspace(EmodelPS.potentialFunc.Workspace.Parameters)
		Ufit2 = eModel2.calculate('current',FullData.current);
        
        RMSE2 = gof2.rmse; % Root mean squares error
        
        % Plotting
        
        figure
        hold on;
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        plot(FullData.current,Ufit2)
        xlabel("$j$ (A)")
        ylabel("$U$ (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Particleswarm')
        hold off;
        
        %% Test with original parameters
%         Emodelfull2 = Emodelfull.copy;
%         Emodelfull2.replaceInWorkspace('alpha',alpha,'j0',j0,'r',r,'j_lim',j_lim)
        eModel3 = eModel.copy;
        eModel3.replaceInWorkspace('alpha',alpha,'j0',j0,'r',r,'j_lim',j_lim)

        
        Ufit = eModel3.calculate('current',FullData.current);
        
        RMSE = sqrt(mean((eModel3.calculate('current',SynData.current(:,1))-SynData.voltage(:,1)).^2)); % Root mean squares error
        
        % Plotting
        
        figure
        hold on;
        errorbar(SynData.current(:,1),SynData.voltage(:,1),SynData.voltage(:,2),SynData.voltage(:,2),SynData.current(:,2),SynData.current(:,2),'o')
        plot(FullData.current,Ufit)
        xlabel("$j$ (A)")
        ylabel("$U$ (V)")
        legend("Data", "Fit", "Location", "Best")
        title('Optimal parameters')
        hold off;


        
    case 2 % Data from Jülich
        
        weights = 'l';
        
        I0fit = nan(4,3);
        alphafit = nan(4,3);
        rfit = nan(4,3);
        RMSE = nan(2,3);
        Rsqrd = nan(2,3);
        
        load('JulichData.mat')
        
        %% Create electrolyzer model objects
        eModel = electrolyzerModel('type','PEM');
        eModel.setInWorkspace(struct('Variables',struct('T',T,'pCat',pH2,'pAn',pO2)));
        eModel.addFuncs('ocv','act','ohm');
        
        for j = 1:3
            
            [Imeas,sortInd] = sort(JulichUI(j).I);
            I = [Imeas JulichUI(j).Istd(sortInd)];
            U = [JulichUI(j).U(sortInd) JulichUI(j).Ustd(sortInd)];
            
            %% Fit
            
            % Non-linear least squares error
            tic;
            [fitParams1,gof1] = eModel.fitUI(U,I,'method','nllse','weights',weights,'plot',true);
            toc
            
            I0fit(1:2,j) = fitParams1.j0;
            alphafit(1:2,j) = fitParams1.alpha;
            rfit(1:2,j) = fitParams1.r;
            
            % Calculate voltage with calculate-method of func-object
            Ufit1 = eModel.calculate('current',Imeas);
            
            RMSE(1,j) = gof1.rmse; % Root mean squares error
            Rsqrd(1,j) = gof1.rsquare; % R^2 of the fit
            
            % Plotting
            
            figure
            hold on;
            errorbar(I(:,1),U(:,1),U(:,2),U(:,2),I(:,2),I(:,2),'o')
            plot(Imeas,Ufit1)
            xlabel("$I$ (A)")
            ylabel("$U$ (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Non-linear least squares error: Jülich ' num2str(j)])
            hold off;
            
            
            % Particleswarm
            tic;
            [fitParams2,gof2] = eModel.fitUI(U,I,'method','ps','weights',weights,'plot',true);
            toc
            
            
            I0fit(3:4,j) = fitParams2.j0;
            alphafit(3:4,j) = fitParams2.alpha;
            rfit(3:4,j) = fitParams2.r;
            
			% Calculate voltage with calculate-method of func-object
            Ufit2 = eModel.calculate('current',Imeas);
            
            RMSE(2,j) = gof2.rmse; % Root mean squares error
            Rsqrd(2,j) = gof2.rsquare; % R^2 of the fit
            
            % Plotting
            
            figure
            hold on;
            errorbar(I(:,1),U(:,1),U(:,2),U(:,2),I(:,2),I(:,2),'o')
            plot(Imeas,Ufit2)
            xlabel("$I$ (A)")
            ylabel("$U$ (V)")
            legend("Data", "Fit", "Location", "Best")
            title(['Particleswarm: Jülich ' num2str(j)])
            hold off;
            
        end
end

%% Test warning message from changing the temperature after fit
eModel.replaceInWorkspace('T',273.15)

%% Test overpotential plotting
j = 0:.01:2;
eModel.showOverpotentials(j)
eModel2.showOverpotentials()
