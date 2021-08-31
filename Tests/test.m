current = 0.1:0.1:2.4;

for i = 1:4
    
    switch i
        case 1
            % Weigh beginning and end of the measured current spectrum
            x = current-mean(current);
        case 2
            % Weigh the beginning of the measured current spectrum
            x = current-min(current);
        case 3
            % Weigh the end of the measured current spectrum
            x = current-max(current);
        case 4
            % Don't apply wieghts
            x = ones(size(current));
    end
    y = x/max(abs(x))*pi/4;
    z = tan(y).^2 + 1;
    
    figure
    plot(z)
end