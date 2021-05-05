function [merged_struct] = mergeStructs(struct_a,struct_b,varargin)

defaultWarnDuplicates = true;

parser = inputParser;
addRequired(parser,'struct_a',@(x) isstruct(x))
addRequired(parser,'struct_b',@(x) isstruct(x))
addParameter(parser,'warn_duplicates',defaultWarnDuplicates,@(x) islogical(x))

parse(parser,struct_a,struct_b,varargin{:});

warn_duplicates = parser.Results.warn_duplicates;

%%
% If one of the structres is empty do not merge
if isempty(struct_a)
    merged_struct=struct_b;
    return
end
if isempty(struct_b)
    merged_struct=struct_a;
    return
end
%%insert struct a
merged_struct=struct_a;

%%insert struct b
overwritten_fields = false;
f_a = fieldnames(struct_a);
f_b = fieldnames(struct_b);
for i = 1:length(f_b)
    if isempty(f_a)
        continue;
    elseif any(f_b{i} == [f_a{:}]) && ~isempty(struct_a.(f_b{i}))
        if isempty(struct_b.(f_b{i}))
            merged_struct.(f_b{i}) = struct_a.(f_b{i});
            continue;
        elseif isstruct(struct_b.(f_b{i})) && isstruct(struct_a.(f_b{i}))
            merged_struct.(f_b{i}) = mergeStructs(struct_a.(f_b{i}),struct_b.(f_b{i}));
        else
            merged_struct.(f_b{i}) = struct_b.(f_b{i});
            overwritten_fields = true;
        end
    else
        merged_struct.(f_b{i}) = struct_b.(f_b{i});
    end
end

if overwritten_fields && warn_duplicates
    warning('Merged structures contained duplicated fields. Values have been overwritten.')
end
end