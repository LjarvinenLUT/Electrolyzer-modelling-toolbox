function MergedStruct = mergeStructs(StructA,StructB,varargin)
% MERGESTRUCTS Merges two structures recursively.
%
%   MergedStruct = MERGESTRUCTS(StructA,StructB)
%   The merged structure contains all the fields from both the original
%   structures. 
%
%   If the same field is contained in both the original structures, the one
%   with a value assigned to it is preferred. If both the structures have a
%   different value assigned to a shared field, then a warning of data loss
%   is issued. This warning can be supressed with a boolean parameter
%   'warn_duplicates' inputed as a name-value pair.
%
%   See also ADDVALUESTOSTRUCT, ISCOMPLETESTRUCT

defaultWarnDuplicates = true;

parser = inputParser;
addRequired(parser,'struct_a',@(x) isstruct(x))
addRequired(parser,'struct_b',@(x) isstruct(x))
addParameter(parser,'warn_duplicates',defaultWarnDuplicates,@(x) islogical(x))

parse(parser,StructA,StructB,varargin{:});

warnDuplicates = parser.Results.warn_duplicates;

%%
% If one of the structres is empty, use the other one
if isempty(fieldnames(StructA))||isempty(StructA)
    MergedStruct=StructB;
    return
end
if isempty(fieldnames(StructB))||isempty(StructB)
    MergedStruct=StructA;
    return
end

%%insert struct a
MergedStruct=StructA;

%%insert struct b
overwrittenFields = false;
overwrittenFieldNames = {};
f_a = fieldnames(StructA);
f_b = fieldnames(StructB);
for i = 1:length(f_b)
    if ismember(f_b{i},f_a) && ~isempty(StructA.(f_b{i}))
        if isempty(StructB.(f_b{i}))
            MergedStruct.(f_b{i}) = StructA.(f_b{i});
            continue;
        elseif isstruct(StructB.(f_b{i})) && isstruct(StructA.(f_b{i}))
            MergedStruct.(f_b{i}) = mergeStructs(StructA.(f_b{i}),StructB.(f_b{i}));
        else
            if ~isequal(StructB.(f_b{i}),StructA.(f_b{i}))
                overwrittenFields = true;
                overwrittenFieldNames = [overwrittenFieldNames,f_b{i}];
            end
            MergedStruct.(f_b{i}) = StructB.(f_b{i});
        end
    else
        MergedStruct.(f_b{i}) = StructB.(f_b{i});
    end
end

if overwrittenFields && warnDuplicates
    fields = overwrittenFieldNames{1};
    for i = 2:length(overwrittenFieldNames)
        fields = [fields ', ' overwrittenFieldNames{i}];
    end
    warningmsg = "Merged structures contained duplicated fields:\n" ...
        + fields ...
        + "\nTheir values have been overwritten by values in the "...
        + "second input structure.";
    warning("mergeStructs:overwritten",warningmsg)
end
end