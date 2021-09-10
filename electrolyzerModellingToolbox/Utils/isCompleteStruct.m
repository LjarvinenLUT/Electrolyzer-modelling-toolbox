function structureComplete = isCompleteStruct(Structure)
% ISCOMPLETESTRUCT Determine recursively if all the fields and subfields in
%   a structure have a value assigned to them.
%
%   See also MERGESTRUCTS, ADDVALUESTOSTRUCT

structureComplete = true; % Set default

fields = fieldnames(Structure);
for i = 1:length(fields)
    if isstruct(Structure.(fields{i}))
        structureComplete = isCompleteStruct(Structure.(fields{i}));
    else
        structureComplete = ~isempty(Structure.(fields{i}));
    end
    
    if ~structureComplete
        break;
    end
end
end