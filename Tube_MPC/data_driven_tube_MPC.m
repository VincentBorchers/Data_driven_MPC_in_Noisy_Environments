clear;
% Load tightened constraints
load('tightened_constraints.mat');


% System definition--------------------------------------------------------

A = [1 1; 0 1];
B = [0; 1];
C = eye(2);
D = [0; 0];

n = 2;
m = 1;
p = 2;

% check for stability
if ~isequal((abs(eig(A)) < ones(n, 1)), ones(n, 1))
    warning('Open-loop system is not stable');
end

% check that system is controllable
if rank(ctrb(A, B)) ~= n
    warning('System is not controllable');
end

% check that system is a minimal representation
sys = ss(A, B, C, D, -1);
sysr = minreal(sys);
if order(sys) ~= order(sysr)
   warning('Order of system and order of minimal realisation differ'); 
end

% Find lag of system
l = 1; % initialises lag
O = []; % initialises obervability matrix / improve by preallocating memory
for i = 1:n
    l = i;
    O = [O; C*A^(i-1)];
    if rank(O) == n
        break
    end
end



% DeePC algorithm ---------------------------------------------------------

% Set algorithm parameters ------------------------------------------------
T_ini = 2;
N = 7;
T = 25;
iterations = 35;
s = 1;

% Define control reference signal 
ref = zeros(N, p);

% Define output cost matrix and control cost matrix
Q = eye(p);
R = eye(m);

% check that T_ini is large enough
if T_ini < l
    warning('T_ini is not large enough');
end

% check that T is large enough
if T < (m+1)*(T_ini+N+n)-1
    warning('T is not large enough');
end


% randomly generate persistently exciting input sequence
u_d = normrnd(0, 1, m, T);
t = 0:1:T-1;

% check that input sequence is indeed persistently exciting of order
% T_ini+N+n
data_first_col = u_d(:, 1:T_ini+N+n);
first_col = data_first_col(:);
last_rows = u_d(:, T_ini+N+n:T);
hankel_u_d = construct_hankel_matrix(first_col, last_rows);
if rank(hankel_u_d) ~= m*(T_ini+N+n)
    warning(['Offline input sequence is not persistently exciting of order T_ini+N+n \n' ...
             'The rank of the Hankel matrix is only %d'], rank(hankel_u_d));
end


% compute offline output sequence
y_d = lsim(sys, u_d.', [], zeros(n,1));

% Build input data Hankel matrix
data_first_col = u_d(:, 1:T_ini+N);
first_col = data_first_col(:);
last_rows = u_d(:, T_ini+N:T);
hankel_u_d_offline = construct_hankel_matrix(first_col, last_rows);

U_p = hankel_u_d_offline(1:m*T_ini, :);
U_f = hankel_u_d_offline((m*T_ini)+1:m*(T_ini+N), :);

% Build output data Hankel matrix
data_first_col = y_d(1:T_ini+N, :).';
first_col = data_first_col(:);
last_rows = y_d(T_ini+N:T, :).';
hankel_y_d_offline = construct_hankel_matrix(first_col, last_rows);

Y_p = hankel_y_d_offline(1:p*T_ini, :);
Y_f = hankel_y_d_offline((p*T_ini)+1:p*(T_ini+N), :);


% Define an arbitrary input sequence that brings the system into a state
% from where the DeePC loop begins.
%u = [0.6 -0.5 0.1];
%u = [0.15 -0.125 0.025];
u = [0 0.4 -1];
y = lsim(sys, u.', [], zeros(n,1));
y = y.';

u_ini_data = u(:, end-T_ini+1:end);
u_ini = u_ini_data(:);

y_ini_data = y(:, end-T_ini+1:end);
y_ini = y_ini_data(:);

% Set lambda for one-norm regularisation of g
lambda_g = 0;

% Build constraint matrix for x-inequality-constraints
A_c_big = zeros(N*p*2, N*p);
for i = 1:N
    A_c_big(1+(i-1)*p*2:1+(i-1)*p*2+2*p-1, 1+(i-1)*p:1+(i-1)*p+p-1) = A_c;
end

% Define parts of the quadratic programming problem that don't change from one
% iteration of the DeePC algorithm to the next one.
A = [eye(T-T_ini-N+1) zeros(T-T_ini-N+1, N*(m+p)) -eye(T-T_ini-N+1) ;...
     -eye(T-T_ini-N+1) zeros(T-T_ini-N+1, N*(m+p)) -eye(T-T_ini-N+1) ; ...
     zeros(N*m,T-T_ini-N+1) eye(N*m) zeros(N*m,N*p+T-T_ini-N+1); ...
     zeros(N*m,T-T_ini-N+1) -eye(N*m) zeros(N*m,N*p+T-T_ini-N+1); ...
     zeros(2*N*p,T-T_ini-N+1+N*m) A_c_big zeros(2*N*p,T-T_ini-N+1)];
Aeq = [[U_p; Y_p; U_f; Y_f] [zeros((m+p)*T_ini, (m+p)*N); -eye((m+p)*N)]];
Aeq = [Aeq zeros((m+p)*(T_ini+N),T-T_ini-N+1)];
H = build_H_matrix(Q, R, N, p, m, T, T_ini); % Build quadratic objective term
H = blkdiag(H, zeros(T-T_ini-N+1, T-T_ini-N+1));
f = build_f_vector(T, T_ini, N, m, Q, ref.'); % Build linear objective term
f = [f; lambda_g*ones(T-T_ini-N+1,1)];

% Build inequality constraint vector b
b = [zeros(2*(T-T_ini-N+1),1); repmat(u_tight_max,N,1);...
     repmat(-u_tight_min,N,1);...
     repmat(b_c,N,1)];


for i = 1:iterations 
    % Define parts of the quadratic programming problem that do change from one
    % iteration of the DeePC algorithm to the next one.
    beq = [u_ini; y_ini; zeros((m+p)*N, 1)];

    % Solve quadratic programming problem
    x = quadprog(H, f, A, b, Aeq, beq);
    disp(i);

    % Extract solution from the x-vector
    g = x(1:T-T_ini-N+1);
    norm(g,1);
    u_fut = reshape(x(T-T_ini-N+1+1:T-T_ini-N+1+m*N), [m,N]);
    y_fut = reshape(x(T-T_ini-N+1+m*N+1:T-T_ini-N+1+m*N+p*N), [p,N]);

    % Update input and output sequences
    u = [u u_fut(:,1:s)]; % take only the first s inputs out of the N-length input sequence
    y = lsim(sys, u.', [], zeros(n,1));
    y = y.';

    % Update u_ini and y_ini
    u_ini_data = u(:, end-T_ini+1:end);
    u_ini = u_ini_data(:);
    y_ini_data = y(:, end-T_ini+1:end);
    y_ini = y_ini_data(:);
end

% Trajectory of nominal system
figure();
scatter(0:1:size(y(:,4:end),2)-1, y(1,4:end));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\bar{x_{1}}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
scatter(0:1:size(y(:,4:end),2)-1, y(2,4:end));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\bar{x_{2}}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
plot([0 0], [1 -1], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot([1 -1], [0 0], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot(y(1,4:end), y(2,4:end), '-o', 'Color', [0 0.4470 0.7410]);
hold on
plot(y(1,4), y(2,4), '.', 'Color', [0 0.4470 0.7410], 'MarkerSize', 30);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$\bar{x}_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$\bar{x}_{2}$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);
xlim([-0.5 0.5]);
ylim([-0.7 0.3]);


% Simulate uncertain system (with additive bounded disturbance)
A = [1 1; 0 1];
K = [-0.4 -1.2];
w_bound = 0.1;
[x_uncertain, u_uncertain] = simulate_uncertain_sys(A, B, u(:,4:end), y(:,4:end), n, m, K, w_bound);

figure();
scatter(0:1:size(x_uncertain,2)-1, x_uncertain(1,:));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$x_{1}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
scatter(0:1:size(x_uncertain,2)-1, x_uncertain(2,:));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$x_{2}$', 'Interpreter', 'latex', 'FontSize', 16, 'Rotation', 0);
grid on;

figure();
plot([0 0], [1 -1], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot([1 -1], [0 0], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot(x_uncertain(1,:), x_uncertain(2,:), '-o', 'Color', [0 0.4470 0.7410]);
hold on
plot(x_uncertain(1,1), x_uncertain(2,1), '.', 'Color', [0 0.4470 0.7410], 'MarkerSize', 30);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$x_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$x_{2}$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);
xlim([-0.5 0.5]);
ylim([-0.7 0.3]);



% Function definitions ----------------------------------------------------

function H = construct_hankel_matrix(first_col, last_rows)
    m = size(last_rows, 1);
    if ~isequal(first_col(end-m+1:end, 1), last_rows(:, 1))
        warning(['Last %d elements of input column do not match first ' ... 
            'column of input rows.\n         Column wins conflict.'], m);
    end
    % preallocate a matrix of the desired size with zeros
    H = zeros(size(first_col, 1), size(last_rows, 2));
    % make last_rows the last rows of H
    H(end-m+1:end, :) = last_rows;
    % make first_col the first column of H
    H(:, 1) = first_col;
    % construct the entire hankel matrix from the first column and the last
    % set of rows
    for i = 2:size(H, 2)
        H(1:end-m, i) = H(m+1:end, i-1);
    end
end

function H = build_H_matrix(Q, R, N, p, m, T, T_ini)
    Q_d= zeros(N*p, N*p);
    R_d = zeros(N*m, N*m);
    for i = 1:N
       Q_d((i-1)*p+1:(i-1)*p+p, (i-1)*p+1:(i-1)*p+p) = Q;
       R_d((i-1)*m+1:(i-1)*m+m, (i-1)*m+1:(i-1)*m+m) = R;
    end
    H = 2*blkdiag(zeros(T-T_ini-N+1, T-T_ini-N+1), R_d, Q_d);
end

function f = build_f_vector(T, T_ini, N, m, Q, r)
    f = zeros(T-T_ini-N+1+N*m, 1);
    for i = 1:N
        f = [f; -2*Q.'*r(:,i)];
    end    
end

function [x, u] = simulate_uncertain_sys(A, B, u_n, x_n, n, m, K, w_bound)
    x = [x_n(:,1) zeros(n, size(u_n,2)-1)];
    u = zeros(m, size(u_n,2));
    for i = 2:size(u_n,2)
        u(:,i-1) = u_n(:,i-1) + K*(x(:,i-1)-x_n(:,i-1));
        x(:,i) = A*x(:,i-1) + B*u(:,i-1) + unifrnd(-w_bound, w_bound, n, 1);
    end
end