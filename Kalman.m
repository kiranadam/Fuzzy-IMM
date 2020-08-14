%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman Filter as per Dr.Wolfgang Koch lecture    %
% Author Kirakumar V. Adam                         %         
%                                                  %
%                                                  %
% xk state vector at time k which is x(k|k)        %
% Pk Covarinace matrix at time k which is P(k|k)   %
% xk_1 state vector such as x(k-1|k-1)             %
% Pk_1 Covariance such as P(k-1|k-1)               %
% F  Evolution matrix at time k-1                  %
% H  Obersvality/Measurement matrix at time k      %
% Q  Process white noise at time k-1               %
% R  Measurement white noise at time k             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[xkk_1,Pkk_1,xk,Pk,nuk,S] = Kalman(xk_1,Pk_1,Zk,F,H,Q,R)

	% Predition 
	xkk_1 = F*xk_1 ;    % state prediction (mean) x(k|k-1)
	Pkk_1 = F*Pk_1*F' + Q ; % covarinace prediction P(k|k-1)

	% Correction/Filtering
	nuk = Zk - H*xkk_1 ; % innovation measurement 
	S   = H*Pkk_1*H' + R; % innovation covarince
	W   = Pkk_1*H'*inv(S); % Kalman Gain

	xk = xkk_1 + W*nuk ; % state filtered (mean) x(k|k)
	Pk = Pkk_1 - W*S*W'; % covariance filtered P(k|k)

return
