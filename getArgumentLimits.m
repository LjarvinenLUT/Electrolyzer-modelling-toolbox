function [lower, upper, start] = getArgumentLimits(argumentList, U, I)
%GETARGUMENTLIMITS gets lower and upper limits and also startpoint of each 
%fitting parameters

lower = zeros(1, length(argumentList));
upper = zeros(1, length(argumentList));
start = zeros(1, length(argumentList));

for i = 1:length(argumentList)
    switch argumentList{i}
        case 'j0'
            lower(i) = 1e-10;
            upper(i) = 1;
            start(i) = 0.001;
        case 'alpha'
            lower(i) = 0;
            upper(i) = 1;
            start(i) = 0.5;
        case 'r'
            lower(i) = 0;
            upper(i) = inf;
            start(i) = 1;
        case 'jL'
            lower(i) = max(I);
            upper(i) = inf;
            start(i) = max(I)+1;
        case 'Uerr'
            lower(i) = -inf;
            upper(i) = inf;
            start(i) = 0;
        otherwise
            error('getArgumentLimits.m: No argument limits could be found for ' + argumentList{i});
    end
end
end

