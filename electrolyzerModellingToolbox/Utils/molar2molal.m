function molality = molar2molal(Molarity,T,solute)
% MOLAR2MOLAL Convert concentration value in 
%   molarity (moles of solute/volume of solution in L) to
%   molality (moles of solute/mass of solvent in kg).
%   Uses iterative method for conversion
%
%   Inputs:
%       Molarity -- concentration as molarity
%       T -- temperature in kelvin
%       solute -- chemical formula of the solute
% 
%    Outputs:
%       molality -- concentration in molality
% 
% See also MOLAL2MOLAR, MOLAL2WTFRAC, WTFRAC2MOLAL

nT = numel(T);
nM = numel(Molarity);

if nT == 1 && nM == 1 % Scalar variables
    trialMolality = [0 25 50];
    % Find the correct molarity iteratively
    maxIter = 100;
    for i = 1:maxIter
        trialMolarity = molal2molar(trialMolality,T,solute,false);
        molarityDif = trialMolarity-Molarity;
        if abs(molarityDif(2)) < 1e-5 % Difference smaller than threshold
            break;
        elseif molarityDif(3) < 0 % Trial concentrations too small
            trialMolality = trialMolality + trialMolality(2);
        elseif molarityDif(2) < 0 % True concentration between trials 2 and 3
            trialMolality = [trialMolality(2) (trialMolality(2)+trialMolality(3))/2 trialMolality(3)];
        elseif molarityDif(2) > 0 % True concentration between trials 1 and 2
            trialMolality = [trialMolality(1) (trialMolality(1)+trialMolality(2))/2 trialMolality(2)];
        end
    end
    if i == maxIter
        warning("Electrolyte density calculation did not converge")
    end

    molality = trialMolality(2);

    % Printing error if concentration or temperature too high:
    solutionDensity(T,molality,solute,true);

    % molality = 1/( solutionDensity(T,solute,Molarity,'molar')./Molarity - M*1e-3 );

else % Vector variables
    if nT == 1 && nM > 1 
        MolarityVec = reshape(Molarity,[nM,1]);
        molalityVec = nan(size(MolarityVec));
        for i = 1:nM
            molalityVec(i) = molar2molal(MolarityVec(i),T,solute);
        end
        molality = reshape(molalityVec,size(Molarity));
    elseif nT > 1 && nM == 1
        TVec = reshape(T,[nT,1]);
        molalityVec = nan(size(TVec));
        for i = 1:nT
            molalityVec(i) = molar2molal(Molarity,TVec(i),solute);
        end
        molality = reshape(molalityVec,size(T));
    elseif size(T) == size(Molarity)
        MolarityVec = reshape(Molarity,[nM,1]);
        TVec = reshape(T,[nT,1]);
        molalityVec = nan(size(MolarityVec));
        for i = 1:nT
            molalityVec(i) = molar2molal(MolarityVec(i),TVec(i),solute);
        end
        molality = reshape(molalityVec,size(T));
    else
        error("Input temperature and molarity vectors have incomparable sizes")
    end


end

end