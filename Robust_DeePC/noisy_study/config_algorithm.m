% Set algorithm parameters ------------------------------------------------
T_ini = 4;
N = 6;
T = 80;
%T = 20;
iterations = 20;
s = 1;
tcLength = n+1; % length of terminal constraint

% Define control reference signal 
% (for now assumed to be the same for each time window)
ref = 0.5*ones(N, 1);

% Define output cost matrix and control cost matrix
Q = 1000*eye(p);
R = eye(m);