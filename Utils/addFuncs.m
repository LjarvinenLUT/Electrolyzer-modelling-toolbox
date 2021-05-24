function newFunc = addFuncs(func1,func2)
% ADDFUNCS Add two func objects together combining their workspaces and
%   using addition for combining their function handles.
%
%   See also FUNC, MERGESTRUCTS 

    if ~isa(func1,'func')||~isa(func2,'func')
        error("Functions to be added should both be of type 'func'")
    end
    
    newFuncEquation = ['@(Workspace) ' func1.getEquationBody ' + ' func2.getEquationBody];
    newFuncHandle = str2func(newFuncEquation);

    NewWorkspace = mergeStructs(func1.Workspace,func2.Workspace);
    
    newFunc = func(newFuncHandle,NewWorkspace);
    
end

