%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIMM Demo                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1; % Delta Timestep

% Evolution matrix  F1
F1 = [eye(2) eye(2).*T eye(2).*(T*T/2);
      zeros(2) eye(2) eye(2).*T;
      zeros(2) zeros(2) eye(2)];

% Track 1
X1 = zeros(6,50);
X1(:,1) = [0 0 50*sin(pi/4) 50*cos(pi/4) 0 -12]';

for i=2:50
    X1(:,i)= F1*X1(:,i-1);
end

abs_v = magnitude(X1(3,:),X1(4,:)); % absolute velocities 1

v1 = max(abs_v); % max velocity 1
u1 = magnitude(X1(5,1),X1(6,1)); % acceleration 1

% Evolution matrix F2
F2 = F1;

% Track 2
X2 = zeros(6,50);
X2(:,1) = [X1(1,50) X1(2,50) 0 0 10 10]';

for i=2:50
    X2(:,i) = F2*X2(:,i-1);
end

abs_v = magnitude(X2(3,:),X2(4,:)); % absolute velocities 2

v2 = max(abs_v); % max velocity 2
u2 = magnitude(X2(5,1),X2(6,1)); % acceleration 2

% Evolution matrix F3
F3 = F2;

% Track 3
X3 = zeros(6,50);
X3(:,1) = [X2(1,50) X2(2,50) 30*sin(pi/3) 30*cos(pi/3) 0 -17]';

for i=2:50
    X3(:,i) = F3*X3(:,i-1);
end

abs_v = magnitude(X3(3,:),X3(4,:)); % absolute velocities 3

v3 = max(abs_v); % max velocity 3
u3 = magnitude(X3(5,1),X3(6,1)); % acceleration 3

% Evolution matrix F4
F4 = F3;

% Track 4
X4 = zeros(6,50);
X4(:,1) = [X3(1,50) X3(2,50) 0 0 22 22]';

for i=2:50
    X4(:,i) = F4*X4(:,i-1);
end

abs_v = magnitude(X4(3,:),X4(4,:)); % absolute velocities 4

v4 = max(abs_v); % max velocity 4
u4 = magnitude(X4(5,1),X4(6,1)); % acceleration 4

% Evolution matrix F4
F5 = F4;

% Track 5
X5 = zeros(6,50);
X5(:,1) = [X4(1,50) X4(2,50) 0 0 24 24]';

for i=2:50
    X5(:,i) = F5*X5(:,i-1);
end

abs_v = magnitude(X5(3,:),X5(4,:)); % absolute velocities 5

v5 = max(abs_v); % max velocity 5
u5 = magnitude(X5(5,1),X5(6,1)); % acceleration 5

% Evolution matrix F4
F6 = F5;

% Track 6
X6 = zeros(6,50);
X6(:,1) = [X5(1,50) X5(2,50) 0 0 26 26]';

for i=2:50
    X6(:,i) = F6*X6(:,i-1);
end

abs_v = magnitude(X6(3,:),X6(4,:)); % absolute velocities 5

v6 = max(abs_v); % max velocity 5
u6 = magnitude(X6(5,1),X6(6,1)); % acceleration 5

X = [X1 X2 X3 X4 X5 X6]; % total track

sigma = 1.8104;

u = [u1 u2 u3 u4 u5 u6]';

mu1 = gaussian(u,u1,sigma); % calculating gaussian  
mu2 = gaussian(u,u2,sigma);
mu3 = gaussian(u,u3,sigma);
mu4 = gaussian(u,u4,sigma);
mu5 = gaussian(u,u5,sigma);
mu6 = gaussian(u,u6,sigma);

[MuF1,record1]= normalizer(mu1); % we get weights in MuF and records  which are gaussian membership 
[MuF2,record2]= normalizer(mu2);
[MuF3,record3]= normalizer(mu3);
[MuF4,record4]= normalizer(mu4);
[MuF5,record5]= normalizer(mu5);
[MuF6,record6]= normalizer(mu6);

% plant noise 
Q1 = eye(6).*(v1^2/2);
Q2 = eye(6).*(v2^2/2);
Q3 = eye(6).*(v3^2/2);
Q4 = eye(6).*(v4^2/2);
Q5 = eye(6).*(v5^2/2);
Q6 = eye(6).*(v6^2/2);

% measurement noise 
R1 = eye(2).*(v1^2);
R2 = eye(2).*(v2^2);
R3 = eye(2).*(v3^2);
R4 = eye(2).*(v4^2);
R5 = eye(2).*(v5^2);
R6 = eye(2).*(v6^2);

% obeservation matrix
H1 = [ eye(2) zeros(2) zeros(2)];
H2 = H1;
H3 = H2;
H4 = H3;
H5 = H4;
H6 = H5;

n = 300;

Z = zeros (2,n);
X_r = X;

for i=1:n
    if i<=50
        X_r(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q1); 
        Z(:,i) = H1*X(:,i) + gauss_rnd(zeros(size(Z,1),1), R1);
    elseif i>50 && i<=100
        X_r(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q2); 
        Z(:,i) = H2*X(:,i) + gauss_rnd(zeros(size(Z,1),1), R2);
    elseif i>100 && i<=150
        X_r(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q3); 
        Z(:,i) = H3*X(:,i) + gauss_rnd(zeros(size(Z,1),1), R3);
    elseif i>150 && i<=200
        X_r(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q4); 
        Z(:,i) = H4*X(:,i) + gauss_rnd(zeros(size(Z,1),1), R4);
    elseif i>200 && i<=250
        X_r(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q5); 
        Z(:,i) = H5*X(:,i) + gauss_rnd(zeros(size(Z,1),1), R5);
    else
        X_r(:,i) = X(:,i) + gauss_rnd(zeros(size(X,1),1), Q6); 
        Z(:,i) = H6*X(:,i) + gauss_rnd(zeros(size(Z,1),1), R6);
   
    end
end

% Weights to update likelihood
modeProb = [0.2 0.2 0.2 0.2 0.1 0.1]';

Transprob = markov_trans(modeProb);


% means and covariances 
m = zeros(6,n); % array of mean for IMM declared
P = zeros(6,6,n); % array of covariance for IMM declared

m_i = zeros(6,6,n);
P_i = zeros(6,6,6,n);


% Model probabilities likelihood array for IMM
MU = zeros(6,1,n);

% State Transtion matrix for IMM
F(:,:,1) = F1;
F(:,:,2) = F2;
F(:,:,3) = F3;
F(:,:,4) = F4;
F(:,:,5) = F5;
F(:,:,6) = F6;

% Process noise for IMM
Q(:,:,1) = Q1;
Q(:,:,2) = Q2;
Q(:,:,3) = Q3;
Q(:,:,4) = Q4;
Q(:,:,5) = Q5;
Q(:,:,6) = Q6;

% Measurement noise for IMM
R(:,:,1) = R1;
R(:,:,2) = R2;
R(:,:,3) = R3;
R(:,:,4) = R4;
R(:,:,5) = R5;
R(:,:,6) = R6;

% Measurement matrices for IMM
H(:,:,1) = H1;
H(:,:,2) = H2;
H(:,:,3) = H3;
H(:,:,4) = H4;
H(:,:,5) = H5;
H(:,:,6) = H6;

% Initial values for 5 models non fuzzy 
xm(:,1) = X(:,1);
xp(:,:,1) = Q1;
xm(:,2) = xm(:,1);
xp(:,:,2) = Q2;
xm(:,3) = xm(:,1);
xp(:,:,3) = Q3;
xm(:,4) = xm(:,1);
xp(:,:,4) = Q4;
xm(:,5) = xm(:,1);
xp(:,:,5) = Q5;
xm(:,6) = xm(:,1);
xp(:,:,6) = Q6;


% Iterative IMM call
for i=1:n
    [MM,PP,modeProb,xm,xp] = IMM(modeProb,Transprob,Z(:,i),F,H,Q,R,xm,xp);
    m(:,i) = MM;
    P(:,:,i) = PP;
    MU(:,:,i) = modeProb;
    m_i(:,:,i) = xm;
    P_i(:,:,:,i) = xp;
     
end


% means and covariances 
mf = zeros(6,n); % array of mean for fuzzy IMM declared
Pf = zeros(6,6,n); % array of covariance for fuzzy IMM declared

mf_i = zeros(6,3,n);
Pf_i = zeros(6,6,3,n);

% Dynamic model for fuzzy
Ff(:,:,1) = F1;
Ff(:,:,2) = F2;
Ff(:,:,3) = F3;

% Dynamic covariance for fuzzy
Qf1(:,:,1) = Q1;
Qf1(:,:,2) = Q2;
Qf1(:,:,3) = Q3;

Qf2(:,:,1) = Q4;
Qf2(:,:,2) = Q5;
Qf2(:,:,3) = Q6;

% Measurement noise for fuzzy
Rf1(:,:,1) = R1;
Rf1(:,:,2) = R2;
Rf1(:,:,3) = R3;

Rf2(:,:,1) = R4;
Rf2(:,:,2) = R5;
Rf2(:,:,3) = R6;

% Observation matrices for fuzzy
Hf(:,:,1) = H1;
Hf(:,:,2) = H1;
Hf(:,:,3) = H1;

% Initial values for 5 models non fuzzy 
xfm(:,1) = X(:,1);
xfp(:,:,1) = Q1;
xfm(:,2) = xm(:,1);
xfp(:,:,2) = Q2;
xfm(:,3) = xm(:,1);
xfp(:,:,3) = Q3;

Transprobf = markov_trans(MuF1);

% Itervative GFIMM call
for i=1:n
    if i<=50
        [MMf,PPf,MuF1,xfm,xfp] = GFIMM(MuF1,Transprobf,Z(:,i),Ff,Hf,Qf1,Rf1,xfm,xfp);
        mf(:,i) = MMf;
        Pf(:,:,i) = PPf;
        %MU(:,:,i) = modeProb;
        mf_i(:,:,i) = xfm;
        Pf_i(:,:,:,i) = xfp;
    elseif i>50 && i<=100
        [MMf,PPf,MuF2,xfm,xfp] = GFIMM(MuF2,Transprobf,Z(:,i),Ff,Hf,Qf1,Rf1,xfm,xfp);
        mf(:,i) = MMf;
        Pf(:,:,i) = PPf;
        %MU(:,:,i) = modeProb;
        mf_i(:,:,i) = xfm;
        Pf_i(:,:,:,i) = xfp;
    elseif i>100 && i<=150
        [MMf,PPf,MuF3,xfm,xfp] = GFIMM(MuF3,Transprobf,Z(:,i),Ff,Hf,Qf1,Rf1,xfm,xfp);
        mf(:,i) = MMf;
        Pf(:,:,i) = PPf;
        %MU(:,:,i) = modeProb;
        mf_i(:,:,i) = xfm;
        Pf_i(:,:,:,i) = xfp;
    elseif i>150 && i<=200
        [MMf,PPf,MuF4,xfm,xfp] = GFIMM(MuF4,Transprobf,Z(:,i),Ff,Hf,Qf2,Rf2,xfm,xfp);
        mf(:,i) = MMf;
        Pf(:,:,i) = PPf;
        %MU(:,:,i) = modeProb;
        mf_i(:,:,i) = xfm;
        Pf_i(:,:,:,i) = xfp;
    elseif i>200 && i<=250
        [MMf,PPf,MuF5,xfm,xfp] = GFIMM(MuF5,Transprobf,Z(:,i),Ff,Hf,Qf2,Rf2,xfm,xfp);
        mf(:,i) = MMf;
        Pf(:,:,i) = PPf;
        %MU(:,:,i) = modeProb;
        mf_i(:,:,i) = xfm;
        Pf_i(:,:,:,i) = xfp;
    else
        [MMf,PPf,MuF6,xfm,xfp] = GFIMM(MuF6,Transprobf,Z(:,i),Ff,Hf,Qf2,Rf2,xfm,xfp);
        mf(:,i) = MMf;
        Pf(:,:,i) = PPf;
        %MU(:,:,i) = modeProb;
        mf_i(:,:,i) = xfm;
        Pf_i(:,:,:,i) = xfp;
        
    end
        
end

% final plotting 
    figure
    plot(X(1,:),X(2,:),'b-',...    % Original road map  
         m(1,:),m(2,:),'r-',...  % IMM filter
         mf(1,:),mf(2,:),'c-',.... % fuzzy IMM
         Z(1,:),Z(2,:),'g-',.....  % observed measurement
         X(1,1),X(2,1),'ko','MarkerSize',10);  %   
    legend('Ground Truth','IMM filter','GFIMM filter','Measured','Start');    
    xlabel('x(t) in units.');
    ylabel('y(t) in units.');
    title('final');
