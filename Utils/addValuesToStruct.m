%%
function Structure = addValuesToStruct(Structure,names,values)

if ~iscell(names)
    names = num2cell(names);
    values = num2cell(values);
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
