%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Multivariate Gaussian Normal Function   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nuk is Innovation mean in Kalman Filter
% S is Innovation Covariance in Kalman Filter
% likelihood is the normal pdf  



function likelihood = gauss_pdf(nuk,S)
    
   % n = length(nuk); % calculating length 
    
    E = -1/2.*(nuk'*inv(S)*nuk); % 
    
   % likelihood = exp(E)/((2*pi)^(n/2)*det(S)^(1/2));
    likelihood = exp(E) / sqrt(2*pi*det(S));
end