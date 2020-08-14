%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMM filter                                              %
% Estimation with Applications to Tracking and Navigation %
% Theory Algorithms and Software                          %
% - Yaakov Bar-Shalom, X.Rong Li                          %
% Page No. 455 - 457.                                     %                                                      
%                                                         %
% modeProb = weights or likelihood                        %
% Transprob = Transition Probability matrix               %
% Z = sensor measurements                                 %
% F = evolution matrix r models                           %
% H = observation matrix r models                         %
% Q = process noise covariance r models                   %
% R = measurement noise covariance r models               %
% xm = means of r models                                  %
% xp = error covariances of r models                      %
% MM = mean of the model                                  %
% PP = covariance of the model                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[MM,PP,modeProb,xm,xp] = IMM(modeProb,Transprob,Z,F,H,Q,R,xm,xp)

	% getting no. of rows and coloumns for transprob matrix
	[r,c] = size(Transprob); 

	% Mixing probabilities calculation
	c_j = Transprob*modeProb; % (rxr)' * (rx1)   equation 11.6.6-8
	MU_ij = Transprob.*(modeProb*(c_j.^(-1))'); % (rxr).*((rx1)*(rx1)') equation 11.6.6-7 

	% Mixed estimates of means and Covariances 
	temp_x = xm * MU_ij;  % mixing mean equation 11.6.6-9
	xk_1k_1 = xm;  % initialization for mixed means difference
	temp_p = xp; % initialization for mixed covariances

	[m,n,d] = size(temp_p); % getting dimension for mixed covariance

	P = zeros(m*n,d); % defination for the mixed covariance 

	for i = 1:c
		xk_1k_1(:,i) = xm(:,i)- temp_x(:,i); % mixed mean difference 
		temp_p(:,:,i) = xp(:,:,i)+ xk_1k_1(:,i)*xk_1k_1(:,i)'; % mixed covarince without MU_ij multiplication     
	end

	P(:) = reshape(temp_p,[m*n d]);

	P = P * MU_ij;  % equation 11.6.6-10

	Pk_1k_1 = reshape(P,m,n,d);

	temp_MU = zeros(d,1);

	% filtering via kalman filter (Call by fucntion) 
	for i = 1:c
		[xkk_1,Pkk_1,xk,Pk,nuk,S] = Kalman(xk_1k_1(:,i),Pk_1k_1(:,:,i),Z,F(:,:,i),H(:,:,i),Q(:,:,i),R(:,:,i));
		xm(:,i) = xk;
		xp(:,:,i) = Pk;
		mu = gauss_pdf(nuk,S); % calculating likelihood 
		temp_MU(i) = mu; % equation 11.6.6-13
	end

	% Mode Probability 
	for i = 1:c
		modeProb(i)= temp_MU(i).*c_j(i)/(c_j'*temp_MU); % equation 11.6.6-15
	end

	% Output Estimate
	MM = xm*modeProb;  % equation 11.6.6.17

	for i= 1:c
		temp_p(:,:,i)=(xp(:,:,i) + (xm(:,i)-MM)*(xm(:,i)-MM)');
	end

	P(:) = reshape(temp_p,[m*n d]);

	PP = zeros(m,n);
	PP(:) = P * modeProb;  % equation 11.6.6-18

	%Pk_1k_1 = reshape(P,m,n,d);
	
end
