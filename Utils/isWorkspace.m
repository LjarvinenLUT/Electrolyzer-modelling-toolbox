function b = isWorkspace(struct)
% ISWORKSPACE Evaluates if a given structure fulfills the requirements for
%   a FUNC class Workspace. Uses the static ISWORKSPACE method of FUNC
%   class. 

b = func.isWorkspace(struct);

end
