%% Electrolyzer modelling tool: Step-by-step usage example of the UI curve fitting functionality
% Usage of the water electrolysis modelling library is mainly based on the
% functionality of an |electrolyzerModel| object. This object is used to
% store all the information reqired for the model including electrolyzer
% type, measured variables, constants, model parameters and
% modelling equations. Simple-to-use methods for fitting and plotting are
% also incorporated.
%
% A basic usage example of performing UI curve fit to measured voltage and
% current data is presented in this file.

%% Constructing the modelling object
% To construct the modelling object one must specify the electrolysis type
% and, in the case of alkaline electrolysis, the chemical formula of the
% electrolyte. Both of these are given as name-value pairs.
eModel = electrolyzerModel('type','alkaline','electrolyte','NaOH');

%%
% Input options for |electrolyzerModel| constructor are:
%
% <html>
% <table>
% <tr border-bottom=1>
%  <th>type</th>
%  <th>electrolyte</th>
% </tr>
% <tr>
%  <td>'pem'</td>
%  <td>'Polymer membrane' (default)</td>
% </tr>
% <tr>
%  <td rowspan="2">'alkaline'</td>
%  <td>'KOH' (default)</td>
% </tr>
% <tr>
%  <td>'NaOH'</td>
% </tr>
% </table>
% </html>
%
% <latex>
% \begin{tabular}{l|l}
% \bf{type} & \bf{electrolyte}\\ \hline
% 'pem' & 'Polymer membrane' (default)\\
% \multirow{2}{*}{'alkaline'} & 'KOH' (default)\\
% & 'NaOH'
% \end{tabular}
% </latex>

%%
% The fresh |electrolyzerModel| object contains the following properties:
eModel.report


%% Defining system variables
% After model construction the simplest way to proceed is to define some
% system variables, like temperature and pressure, which are needed for
% calculating the Nernst potential for the system. For alkaline
% electrolysis the required variables are:
%
% * Temperature, T (in kelvin)
% * System pressure, ps (in bara)
% * Electrolyte concentration, either as
%   molality (in mol of solute/kg of solvent),
%   Molarity (in mol of solute/L of solution) or
%   wtfrac (in mass of solute/mass of solution)
%
% For PEM, on the other hand, the required variables would be:
%
% * Temperature, T (in kelvin)
% * Anode pressure, pAn (in bara)
% * Cathode pressure, pCat (in bara).
%
T = 273.15 + 70; % Temperature in kelvin
ps = 3; % System pressure in bara

%%
% The alkaline models require concentration as molality, but often weight
% fractions are more commonly known. Therefore the |electrolyzerModel|
% class contains methods for automatic conversion between molality and 
% weight fractions. These conversions use the known molar mass of the given
% electrolyte salt. Weight fraction can be given both as percents (assumed
% when input > 1) or directly as a fraction (input < 1).
wtfrac = 30; % concentration as weight percentage

%%
% Variables have to be provided to the model as a Workspace structure (see
% the documentation of |func| for more information). A Workspace structure
% contains substructures for Variables, Parameters and Constants. For our
% purpose we have to provide the previously defined variables:
Workspace = struct('Variables',struct('T',T,'ps',ps,'wtfrac',wtfrac))

%%
% To check compatibility with the Workspace requirements, one can use the
% static method |func.isWorkspace| of the |func| class.
func.isWorkspace(Workspace)

%%
% The workspace-compatible structure can now be set for the electrolyzer
% model:
eModel.setInWorkspace(Workspace)
eModel.viewWorkspace;

%% Defining the modelling function
% Once the modeling object has been created and variables set, it is
% possible to use helper method |addPotentials| to construct the full
% overpotential function.
eModel.addFuncs('nernst','ohmic','activation')

%%
% The input above provides the default versions of each potential term
% using the information already included in the model object. Alternatively
% one can input non-default potential terms by providing any number of
% |func| objects that are either created by one of the functions
% |nernst|, |ohmic|, |activation| or |concentration|, or custom built
% by the user. Most of them existing potential functions have alternative
% model options available that can be used in the modelling by providing
% the |func| -objects directly (see the documentation of the respective
% functions). One can also provide any other |func| object for including
% some custom potential term that has not been listed in this documentation
% file.
%
% Variable amount of potential terms can be added and the order for them
% does not matter. Some methods like the |fitUI| and |showUI|
% cannot be called before potentials are added as there would be nothing to
% fit or plot in that case.
%
% The |func| objects provided for the system all introduced their necessary
% constants and parameters to the Workspace structure
eModel.report;

%% Copying the modelling object
% If one wants to copy the created electrolyzer model before further
% operations, direct reassigning is not going to do the trick. As the class
% |electrolyzerModel| inherits class |handle|, changes to a reassigned
% object will affect the original one and vice versa. Therefore, the class
% has a method |copy| included that creates a separate object that includes
% all the information from its parent but breaks the link between the child
% and the parent objects.
eModel2 = eModel.copy;
if isequal(eModel,eModel2)
    disp("Copy (eModel2) is the same object as the original (eModel)")
else
    disp("Copy (eModel2) is separate from the original (eModel)")
end


%% Fitting
% To fit the model to existing UI measurement data one can use the method
% |fitUI| of |electrolyzerModel| class.
%
% Some example data can be found from folder TestData:
load('AlkaliData.mat')
voltageData = AlkaliUI.U;
currentData = AlkaliUI.j;
temperatureData = AlkaliUI.T(:,1);
pressureData = AlkaliUI.P(:,1);

%%
% Alternatively one can create synthetic UI curve for testing purpose using
% special function |createSyntheticUI|, which enables user to set the
% amount of data points, defaulted at 20. The user can also induce normally
% distributed measurement error by defining the number of measurements for
% each data point and the measurement error as a fraction of the reading.
% For more detailed explanation of the function, see its
% <matlab:doc('createSyntheticUI') documentation>.

%%
% Let's replace the preset temperature from the electrolyzer model with the
% measured temperature vector
eModel.replaceInWorkspace('T',temperatureData,'ps',pressureData);

%%
% Let's choose particle swarm as our fitting method.
method = "PS";
%%
% Alternatively one could use the Non-Linear Least Squares Error regression
% by calling for |"NLLSE"|.
%
% Weighting of the low current values is enabled with the following option:
weights = "l";
%%
% The weights are added to improve parametrization of the activation
% overpotential, whose effect is post prominent in the lower current
% densities. Now that mass transfer effects are not present in the data to
% be fitted, we do not weigh the higher current densities, which could be
% done by adding letter "h" to the weights call. More in-detail description
% of the options can be found from the <matlab:doc('fitUI') documentation
% of the function |fitUI|>.

[fitParams,gof] = eModel.fitUI(voltageData,currentData,'method',method,'weights',weights);

%%
% *NOTE:* The fitting tool doesn't consider the units of the measured data
% but the user has to keep in mind the used units. Some fitting
% parameters are sensitive to units, for example resistance |r| and
% exchange current density |j0|.

%% Viewing the results
% To see the results, one can use the |showUI| method to perform a quick
% automated plot.
eModel.showUI

%% 
% The parameter values and their uncertainty (standard deviation) can be
% seen from the output of |fitUI|
disp(fitParams)
%%
% or by calling |electrolyzerModel.getParams| method.
disp(eModel.getParams)


%%
% Some goodness of fit values are stored in |fitUI| output |gof|: 
%
% * ssr: Square Sum Residuals
% * rmse: Root Mean Square Error
% * rsquare: the R^2 value of the fit
%
disp(gof)

%% Using fit parameters
% To calculate the voltage values based on the just fitted UI curve, one
% can use the |calculate| method of |electrolyzerModel|.
Ufit = eModel.calculate('current',currentData);

%%
% The variables still missing from the Workspace have to be provided as
% name value pairs. The variables given as input to the |calculate| method
% are preferred over the values with the same name already contained in the
% Workspace. This way the model can be used to calculate cell voltage based
% on the UI curve in different conditions than where the curve used for
% parametrization was measured.

