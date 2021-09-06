%% Creating custom potential term as a |func| object
% If the user cannot find the desired models from the modelling functions
% provided with the toolbox, they can define their own potential terms by
% creating custom |func| objects. To create a |func| object one needs a
% |Workspace| structure and the function handle that uses *only* the
% |Workspace| as its input. The |Workspace| structure should contain all
% the required variables, constants and coefficients (at least empty
% fields) in the |Variables|, |Constants| and |Coefficients| substructures,
% respectively. If there are any dependencies between any of the
% abovementioned parameters, these should be stored as strings in the
% |Dependencies| substructure. The values of the parameters can be left
% empty for later addition (recommended when there are multiple potential
% terms using the same parameters) or they can be preset in the creation
% process.
%
% Let's create a simple custom |func| object as an example. The easiest way
% is to build the function handle first to see which parameters are
% required for the |Workspace|
funcHandle = @(Workspace) Workspace.Constants.a*Workspace.Variables.x.^Workspace.Coefficients.b - Workspace.Variables.t;

%%
% The function handle above represents function
%
% $$a x^b - t.$$
%
% For this potential term the |Workspace| contains one constant $a$, two
% variables $x$ and $t$, and one coefficient $b$.
%
% The substructures needed for the |Workspace| can now be created:
Constants = struct('a',2);
Variables = struct('x',[],'t',[]);
Coefficients = struct('b',[]);
Dependencies = struct('t','Workspace.Variables.t = sin(5*pi.*Workspace.Variables.x);');

%%
% Here theconstant $a$ has been assigned value 2 in the creation process
% but the other parameters are left empty for now. Also the variable $t$ is
% defined to be dependent on the value of variable $x$ through function
% $t = \sin(5\pi x)$.
%
% After the substructures have been defined, the |Workspace| structure can
% be built:
Workspace = struct('Variables',Variables,...
                   'Constants',Constants,...
                   'Coefficients',Coefficients,...
                   'Dependencies',Dependencies);

%%
% The |func| object itself can now be created by calling the constructor
% method:
theFunc = func(funcHandle,Workspace);

%%
% To add missing parameters to the object, the method |func.replaceParams|
% has to be used.
theFunc.replaceParams('x',(0:0.01:1)');
disp(theFunc.viewWorkspace)

%%
% Alternatively the the |func| object could now be included to an
% |electrolyzerModel| object with |electrolyzerModel.addPotentials| method.
% Parameters can be replaced also with the help of
% |electrolyzerModel.replaceParams| method. 
%
% For determining the value for parameter $b$ through curve fitting, its
% value has to be left empty. For this example, we assume the fit has been
% performed and the result has been:
b = 3.25;
theFunc.replaceParams('b',b);

%%
% With all the parameters defined, the value of the function can be
% calculated using method |calculate|:
y = theFunc.calculate();

figure
plot(theFunc.Workspace.Variables.x,y)
xlabel('x')
ylabel('y')

%% Reason for using |func| objects
% This example doesn't demonstrate the reason behind the usage |func|
% objects instead of directly using the function handles for fitting.
% The special features of |func| class enable three main aspects:
%
% # It enables the addition of two |func| objects with a static method
% |func.add| so that the function handle is maintained in readable form.
% # Adding two func objects unites the Workspaces of these two enabling
% later modifications of the parameters.
% # Using only one structure as the parameter for the function handle
% removes the need for correct input order for the function handle, as the
% function handle calls only the parameters it needs and uses them only in
% the intended way.
%
% The flexibility of the Electrolyzer Modelling Toolbox is mainly due to
% the ability of combining multiple potential functions in desired way
% still maintaining the ability to perform parametrization through curve
% fitting. At the same time, the function stays fully readable and can be
% modified as needed.
