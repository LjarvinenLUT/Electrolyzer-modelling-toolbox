function structureComplete = isCompleteStructure(Structure)
% Check recursively if all the fields and subfields in a 
% structure have a value.

structureComplete = true; % Set default

fields = fieldnames(Structure);
for i = 1:length(fields)
    if isstruct(Structure.(fields{i}))
        structureComplete = isCompleteStructure(Structure.(fields{i}));
    else
        structureComplete = ~isempty(Structure.(fields{i}));
    end
    
    if ~structureComplete
        break;
    end
end
end