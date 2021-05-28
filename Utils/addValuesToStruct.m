function Structure = addValuesToStruct(Structure,varargin)
% ADDVALUESTOSTRUCT Adds user specified values to the user specified fields
%   in the given structure. The function doesn't add fields but only adds
%   values to existing ones.
%
%   See also MERGESTRUCTS, ISCOMPLETESTRUCT

if isempty(varargin{1})
    return;
elseif length(varargin) == 1 && isstruct(varargin{1}) % structure input
    names = fieldnames(varargin{1});
    values = struct2cell(varargin{1});
elseif length(varargin) == 2 && iscell(varargin{1}) % array inputs
    names = varargin{1};
    if  isnumeric(varargin{2})
        values = num2cell(varargin{2});
    else
        values = varargin{2};
    end
elseif mod(nargin,2) % name-value pair inputs
    names = {};
    values = {};
    for i = 1:length(varargin)/2
        names{i} = varargin{2*i-1};
        values{i} = varargin{2*i};
    end
else
    error('Input is not valid. Provide the arguments in either a single structure, a cell array of names and an array of values, or as name-value pairs.')
end

if any(isstruct(values))
    error('addValueToStructure is only for adding single field values recursively in a structure')
else
    fn = fieldnames(Structure);
    for i = 1:length(fn)
        if isstruct(Structure.(fn{i}))
            Structure.(fn{i}) = addValuesToStruct(Structure.(fn{i}),names,values);
        elseif ismember(fn{i},names)
            Structure.(fn{i}) = values{strcmp(fn{i},names)};
        end
    end
    
end
end
