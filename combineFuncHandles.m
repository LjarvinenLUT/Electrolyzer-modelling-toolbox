function new_func_handle = combineFuncHandles(func_handles)

z = 0;

for i = 1:length(func_handles)
    % Parse argument names from function handles to cell arrays
    if isa(func_handles{i},'function_handle')
        args = getFunctionArguments(func_handles{i},'omitCurrent',false);
        symargs = num2cell(sym(args));
        zi = func_handles{i}(symargs{:});
    else
        zi = func_handles{i};
    end
    % Combined symbolic result
    z = z+zi;
end

temp_func_handle = matlabFunction(z);
args = getFunctionArguments(temp_func_handle,'omitCurrent',false);
Current_pos = find(strcmp(args,'Current'));
sortindex = [1:Current_pos-1 Current_pos+1:length(args) Current_pos];
args = args(sortindex);
new_func_handle = matlabFunction(z,'Vars',args);
end