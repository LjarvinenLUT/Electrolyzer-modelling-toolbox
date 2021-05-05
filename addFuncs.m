function [new_func] = addFuncs(func_a,func_b)
% Merge two func objects to a new one containing all the information from
% both the original func objects
    new_func_name = [combineHeader(func_a,func_b) body(func_a) '+' body(func_b)];
    new_workspace = mergeStructs(func_a.workspace,func_b.workspace);
    new_func = func(str2fun(new_func_name),new_workspace);
end

function header = combineHeader(func_a,func_b)
    
end

function body = body(func)

end