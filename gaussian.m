function[mu] = gaussian(x,c,sigma)
 
    n = length(x);
    mu = zeros(n,1);
    
    for i =1:n
        e = -(x(i)-c)^2/(2*sigma^2);
        mu(i) = exp(e);
    end
    
end