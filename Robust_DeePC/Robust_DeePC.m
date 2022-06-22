clear;

% System definition--------------------------------------------------------

    
% Set system parameters ---------------------------------------------------
mass = 1;
spring_k = 1;
damping_c = 3;

% Define a continuous-time LTI system in state-space form 
% (mass-damper-spring system)
zeta = damping_c/(2*sqrt(spring_k*mass));
A_c = [0 1; -spring_k/mass -damping_c/mass];
B_c = [0; 1/mass];
C_c = [1 0];
D_c = [0];

sysc = ss(A_c, B_c, C_c, D_c);

% convert to discrete-time system using zero-order hold and Ts=0.1
sysd = c2d(sysc, 0.1);

[A, B, C, D] = ssdata(sysd);

n = size(A, 2);
m = size(B, 2);
p = size(C, 1);

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
T_ini = 4;
N = 7;
T = 80;
iterations = 36;
s = 1;
tcLength = n+1; % length of terminal constraint

% Define control reference signal 
% (for now assumed to be the same for each time window)
ref = 0.5*ones(N, 1);
x_0 = [0; 0];

% Define output cost matrix and control cost matrix
Q = 1000*eye(p);
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
y_d_power = sum(abs(y_d).^2, 'all')/size(y_d,1);

% Add Gaussian noise
%noise_std = 0.004;%sqrt(0.01*y_d_power);
noise_std = 0.012;
y_d = y_d + normrnd(0, noise_std, T, 1);


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
u = [1 1.5 2 1];
y = lsim(sys, u.', [], x_0);
y = y.';

u_ini_data = u(:, end-T_ini+1:end);
u_ini = u_ini_data(:);

y_ini_data = y(:, end-T_ini+1:end);
y_ini = y_ini_data(:);

% Set lambda for one-norm regularisation of g
%lambda_g = 100;
lambda_g = 300;

% Define parts of the quadratic programming problem that don't change from one
% iteration of the DeePC algorithm to the next one.
A = [eye(T-T_ini-N+1) zeros(T-T_ini-N+1, N*(m+p)) -eye(T-T_ini-N+1) ;...
     -eye(T-T_ini-N+1) zeros(T-T_ini-N+1, N*(m+p)) -eye(T-T_ini-N+1)];
b = zeros(2*(T-T_ini-N+1),1);
Aeq = [[U_p; Y_p; U_f; Y_f] [zeros((m+p)*T_ini, (m+p)*N); -eye((m+p)*N)]];
Aeq = [Aeq ; [zeros(tcLength, T-T_ini-N+1+(m+p)*N-tcLength) eye(tcLength)]];
Aeq = [Aeq zeros((m+p)*(T_ini+N)+tcLength,T-T_ini-N+1)];
H = build_H_matrix(Q, R, N, p, m, T, T_ini); % Build quadratic objective term
H = blkdiag(H, zeros(T-T_ini-N+1, T-T_ini-N+1));
f = build_f_vector(T, T_ini, N, m, Q, ref.'); % Build linear objective term
f = [f; lambda_g*ones(T-T_ini-N+1,1)];


for i = 1:iterations 
    % Define parts of the quadratic programming problem that do change from one
    % iteration of the DeePC algorithm to the next one.
    beq = [u_ini; y_ini; zeros((m+p)*N, 1); ones(tcLength,1)*ref(1,1)];

    % Solve quadratic programming problem
    % options = optimoptions('quadprog', 'MaxIterations', 200);
    % x = quadprog(H, f, A, b, Aeq, beq,[],[],[],options);
    x = quadprog(H, f, A, b, Aeq, beq);
    disp(i);

    % Extract solution from the x-vector
    g = x(1:T-T_ini-N+1);
    norm(g,1);
    u_fut = reshape(x(T-T_ini-N+1+1:T-T_ini-N+1+m*N), [m,N]);
    y_fut = reshape(x(T-T_ini-N+1+m*N+1:T-T_ini-N+1+m*N+p*N), [p,N]);

    % Update input and output sequences
    u = [u u_fut(:,1:s)]; % take only the first s inputs out of the N-length input sequence
    y = lsim(sys, u.', [], x_0);
    y = y.';

    % Update u_ini and y_ini
    u_ini_data = u(:, end-T_ini+1:end);
    u_ini = u_ini_data(:);
    y_ini_data = y(:, end-T_ini+1:end);
    y_ini = y_ini_data(:);
end

error = sum(abs(y(:,T_ini+1:end)-ref(1,1)*ones(1, size(y(:,T_ini+1:end), 2))), 'all');
disp(error);


% Plot system output obtained from DeePC
figure();
scatter(0:1:size(y,2)-1, y);
hold on
plot(0:1:size(y,2)-1, ref(1,1)*ones(1,size(y,2)),'LineStyle', '--', 'Color', [0.545 0 0]);
ax = gca;
ax.FontSize = 16;
xlabel('\textbf{time step}', 'Interpreter','latex','FontSize',20);
ylabel('\boldmath$y$', 'Interpreter', 'latex', 'FontSize',20, 'Rotation',-360);
legend('y', 'reference', 'Location', 'east', 'Fontsize', 20, 'Interpreter', 'latex');
grid on;

% Plot system input obtained from DeePC
figure();
scatter(0:1:size(u,2)-1, u);
ax = gca;
ax.FontSize = 16;
xlabel('\textbf{time step}', 'Interpreter','latex','FontSize',20);
ylabel('\boldmath$u$', 'Interpreter', 'latex', 'FontSize',20, 'Rotation',-360);
legend('u', 'Location', 'east', 'Fontsize', 20, 'Interpreter', 'latex');
grid on;






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