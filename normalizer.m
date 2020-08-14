%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalizer as eqn 14 IEE Proc-Radar, Sonar Navig.. Vol. 145, No. 6, %
% December 1998 Irwin                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [MuF,record] = normalizer(mu)

    n = length(mu);
    
    num_dig = 2;
    j = 1 ;
    
    for i=1:n
        temp = mu(i)/sum(mu);
        temp = round(temp*(10^num_dig))/(10^num_dig);
        
        if temp > 0.0009
            MuF(j,1) = temp;
            record(j,1) = i;
            j=j+1;
        end
    end
            
end