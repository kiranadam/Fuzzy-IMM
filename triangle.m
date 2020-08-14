function [mu] = triangle(x,c,w)

    if abs(x-c)>=w/2
        mu = 0;
    elseif abs(x-c)<w/2
        mu = 1- abs(x-c)/(w/2);
    end

end