%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMM Demo                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	T = 1; % timestep 

	F1 = [eye(2) eye(2).*T eye(2).*(T^2/2);
		  zeros(2) eye(2) eye(2).*T;
          zeros(2) zeros(2) eye(2)];


	X1 = zeros(6,30);
	X1(:,1) = [0 500 5*cos(pi/3) 5*sin(pi/3) 0 -9.81]';

	for i=2:30
		X1(:,i) = F1*X1(:,i-1);
	end

	F2 = F1;

	X2 = zeros(6,20);
	X2(:,1) = [X1(1,30) X1(2,30) 0 0 15 15]';

	for i=2:20
		X2(:,i) = F2*X2(:,i-1);
	end

	X = [X1 X2];

	plot(X(1,:),X(2,:),'-b');

	n = 50;

	Q1 = eye(6).*((285)^2);
	Q2 = eye(6).*(280^2);
  
	R1 = eye(2).*((285)^2);
	R2 = eye(2).*((280)^2);

	H1 = [eye(2) zeros(2) zeros(2)];
	H2 = H1;

	Z = zeros(2,n);
	X1 = X;

	for i=1:n
		if i<=30
			X1(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q1); 
			Z(:,i) = H1*X1(:,i) + gauss_rnd(zeros(size(Z,1),1), R1);
		else
			X1(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q2); 
			Z(:,i) = H2*X1(:,i) + gauss_rnd(zeros(size(Z,1),1), R2);
		end
	end

	% transition probability
	Transprob = [ 0.8  0.2;
                  0.1  0.9];

	% Weights to update likelihood
	modeProb = [0.5 0.5]';

	% means and covariances 
	m = zeros(6,n); % array of mean for IMM declared
	P = zeros(6,6,n); % array of covariance for IMM declared

	m_i = zeros(6,2,n);
	P_i = zeros(6,6,2,n);

	% Model probabilities likelihood array for IMM
	MU = zeros(2,1,n);

	% State Transtion matrix for IMM
	F(:,:,1) = F1;
	F(:,:,2) = F2;

	% Process noise for IMM
	Q(:,:,1) = Q1;
	Q(:,:,2) = Q2;

	% Measurement noise for IMM
	R(:,:,1) = R1;
	R(:,:,2) = R2;

	% Measurement matrices for IMM
	H(:,:,1) = H1;
	H(:,:,2) = H2;

	% Initial values for 2 models
	xm(:,1) = X(:,1);
	xp(:,:,1) = Q1;
	xm(:,2) = xm(:,1);
	xp(:,:,2) = Q2;

	% Iterative IMM call
	for i=1:n
		[MM,PP,modeProb,xm,xp] = IMM(modeProb,Transprob,Z(:,i),F,H,Q,R,xm,xp);
		m(:,i) = MM;
		P(:,:,i) = PP;
		MU(:,:,i) = modeProb;
		m_i(:,:,i) = xm;
		P_i(:,:,:,i) = xp;
    
		% plotting 
		figure(2)
		plot(X(1,1:i),X(2,1:i),'b-',...    % Original road map  
			m(1,1:i),m(2,1:i),'r-',...  % IMM filter
			Z(1,1:i),Z(2,1:i),'g-',.....  % observed measurement
			X(1,1),X(2,1),'ko','MarkerSize',10);  %   
		legend('Ground Truth','IMM filter','Measured','Start');    
		xlabel('x(t) in units.');
		ylabel('y(t) in units.');
		title('Simulation');
	end


