%% Electrolyzer modelling tool: Step-by-step usage example of the UI curve fitting functionality
% Usage of the water electrolysis modelling library is mainly based on the
% functionality of an |electrolyzerModel| object. This object is used to
% store all the information reqired for the model including electrolyzer
% type, measured variables, constants, model coefficients and
% modelling equations. Simple-to-use methods for fitting and plotting are
% also incorporated.
%
% A basic usage example of performing UI curve fit to measured voltage and
% current data is presented in this file.

%% Constructing the modelling object
% To construct the modelling object one must specify the electrolysis type
% and, in the case of alkaline electrolysis, the chemical formula of the
% electrolyte. Both of these are given as name-value pairs.
eModel = electrolyzerModel('type','alkaline','electrolyte','KOH');

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
% * Electrolyte molality, m (in mol of solute/kg of solvent).
%
% For PEM, on the other hand, the required variables would be:
%
% * Temperature, T (in kelvin)
% * Anode pressure, pAn (in bara)
% * Cathode pressure, pCat (in bara).
%
T = 273.15 + 80; % Temperature in kelvin
ps = 3; % System pressure in bara

%%
% The alkaline models require concentration as molality, but often weight
% fractions are more commonly known. Therefore the |electrolyzerModel|
% class contains methods for conversion between molality and weight
% fractions. These conversions use the known molar mass of the given
% electrolyte salt. Weight fraction can be given both as percents (assumed
% when input > 1) or directly as a fraction (input < 1).
wtfrac = 30; % concentration as weight percentage
m = eModel.wtfrac2mol(wtfrac)

%%
% Variables have to be provided to the model as a Workspace structure (see
% the documentation of |func| for more information). A Workspace structure
% contains substructures for Variables, Coefficients and Constants. For our
% purpose we have to provide the previously defined variables:
Workspace = struct('Variables',struct('T',T,'ps',ps,'m',m))

%%
% To check compatibility with the Workspace requirements, one can use the
% static method |func.isWorkspace| of the |func| class.
func.isWorkspace(Workspace)

%%
% The workspace-compatible structure can now be set for the electrolyzer
% model:
eModel.setParams(Workspace)
eModel.viewWorkspace;

%% Defining the modelling function
% Once the modeling object has been created and variables set, it is
% possible to use helper method |addPotentials| to construct the full
% overpotential function.
eModel.addPotentials('nernst','ohmic','activation','concentration')

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
% constants and coefficients to the Workspace structure
eModel.report;

%% Copying
% If one wants to copy the created electrolyzer model before further
% operations, direct reassigning is not going to do the trick. As the class
% |electrolyzerModel| inherits class |handle|, changes to a reassigned
% object will affect the original one and vice versa. Therefore, the class
% has a method |copy| included that creates a separate object that includes
% all the information from its parent but breaks the link between the child
% and the parent objects.
eModel2 = eModel;
isequal(eModel,eModel2) 

eModelCopy = eModel.copy;
isequal(eModel,eModelCopy)




%% Fitting
% To fit the model to existing UI measurement data one can use the method
% |fitUI| of |electrolyzerModel| class.
%
% Some example data can be found from folder TestData:
load('AlkaliData.mat')
voltageData = AlkaliUI.U;
currentData = AlkaliUI.j;
temperatureData = AlkaliUI.T(:,1);

%%
% Let's replace the preset temperature from the electrolyzer model with the
% measured vector
eModel.replaceParams('T',temperatureData);

%%
% Let's choose particle swarm as our fitting method and enable weighting of
% the beginning and the end of the dataset. More in-detail description of
% the options can be found from the documentation of the function |fitUI|.
method = "PS";
weights = "default";
[fitParams,gof] = eModel.fitUI(voltageData,currentData,'method',method,'weights',weights)

%%
eModel.showUI
