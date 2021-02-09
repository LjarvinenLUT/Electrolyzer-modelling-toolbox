% 

y1 = @computeSquare;
y2 = @multiply;
% x = 5;
% a = 2;
y = @(x,a)y1(x) + y2(x,a)


function y = computeSquare(x)
    y = x.^2;
end

function y = multiply(x,a)
    y = x.*a;
end
