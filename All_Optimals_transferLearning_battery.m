clc;
clear;
close all;

%% Define general parameters
nStates = 3;
nInputs = 1;
nRealization = 1;
nSample = 100000;


%% Time sequence
T0 = 0;
Ts = 1;
Tf = 1000;

t = T0:Ts:Tf;
t_ = T0:Ts:Tf+Ts;

nSteps = numel(t);
%% Noise generation
mu = [0 0 0];
varNoise = 0.0005;
coVarNoise = 0;
coVarMatrix = varNoise * eye(nStates);
sigma = coVarMatrix;

noise = noise_generator(mu,sigma,nStates,nSample,nRealization);


%% Define the parameters related to the source system

At = [0.0678, 0, 0;
      0, 0.9609, 0
      0, 0, 0.98];

Bt = [0.0257;
      0.2850;
     -0.0250];
%% The difference between the actual system and the target system
ep_AB = 0.08;
delta_A =  ep_AB*ones(nStates);
delta_B = ep_AB*ones(nStates,nInputs);
%% Define the parameters of the target system

As = At - delta_A/2;
Bs = Bt - delta_B/2;

%% Define parameters for uncertain system 

SourceSys = [As, Bs]';
A_us = eye(nStates+nInputs);
B_us = -[As, Bs]';
C_us = SourceSys'*SourceSys - ep_AB^2*eye(nStates);

mat_us = -[C_us, B_us';
           B_us, A_us ];

padded_mat_us = padarray(mat_us, [nStates nStates], 0, 'post');
%% Define parameters for collected data

% first we need to collect some data from the target system
N = nStates + 1;

X = zeros(nStates,N+1);
X(:,1) = [1 1 1]';
U0 = 5*randn(nInputs,N);
W0 = noise(1:N,1:nStates)';

for kk = 1:N

    X(:,kk+1) = At*X(:,kk) + Bt*U0(:,kk) + W0(:,kk);
    % X(:,kk+1) = At*X(:,kk) + Bt*U0(:,kk); %no noise

end

X0 = X(:,1:N);
X1 = X(:,2:N+1);

data = [U0; X0];
data_use = [X0; U0];
% Check rank
if(rank(data)<min(size(data)))

    errordlg('Data is not full row rank!');

end

% Now we define the parameters for the data based uncertainty
A_dt = data_use*data_use';
B_dt = -data_use*X1';
C_dt = -sigma + X1*X1';

mat_dt = -[  C_dt, B_dt';
             B_dt, A_dt ];

padded_mat_dt = padarray(mat_dt, [nStates nStates], 0, 'post'); 


%% Define parameters for optimal ellipsoid that falls at the intersection of the two above ellipsoids

W_x = [1,0,0;
    0,1,0;
    0,0,1];

W_u = 1;


Qo = sdpvar(nInputs,nStates);
Po = sdpvar(nStates,nStates);
L = sdpvar(nStates+nInputs, nStates+nInputs);


A_lm = [Po-eye(nStates), Qo',zeros(nStates);
        Qo, zeros(nInputs),-Qo;
        zeros(nStates),-Qo', Po];
B_lm = zeros(nStates+nInputs, nStates);
C_lm = -Po+eye(nStates);

mat_lm = -[C_lm, zeros(nStates), zeros(nStates,nInputs), zeros(nStates);
         zeros(nStates)',Po, Qo',zeros(nStates);
          zeros(nStates,nInputs)',Qo, zeros(nInputs),-Qo;
         zeros(nStates)',zeros(nStates),-Qo', -Po];


tau_dt = sdpvar(1,1);
tau_us = sdpvar(1,1);
miu_dt =sdpvar(1,1);
miu_us =sdpvar(1,1);
gamma = sdpvar(1,1);


constraintO_1 = - miu_dt * padded_mat_dt -miu_us * padded_mat_us>= 0.0;
constraintO_2 = mat_lm - tau_dt*padded_mat_dt - tau_us*padded_mat_us >=0;
constraintO_3 = tau_us>=0.00000000001;
constraintO_4 = tau_dt>=0.00000000001;
constraintO_5 = L - [0 , Qo;
                    Qo', Po]>=0;
constraintO_6 = trace(W_x*Po) +trace(W_u*L)<=gamma;
constraintO_7 = Po>=eye(nStates);

constraints_o_ =   constraintO_1 + constraintO_2 + constraintO_3+ constraintO_4  + constraintO_5 + constraintO_6 + constraintO_7 ;
cost_o_ = gamma;

ops = sdpsettings('solver','mosek','verbose',0);
diagnostics = optimize(constraints_o_, cost_o_, ops);

if(diagnostics.problem == 0)

    P_opt_o = value(Po);
    Q_opt_o = value(Qo);
    tau_dt_opt = value(tau_dt);
    % lambda_opt = value(lambda);

else

    fprintf('Error! Solution is infeasible!\n\n');

end

K_opt = Q_opt_o*inv(P_opt_o);



x0 = [15 10 -25]';


%% transfer learning lqr

% Main loop
x = zeros(nRealization,nStates,nSteps+1);

for ii = 1:nRealization

    % Simulation
    x(ii,:,1) = x0';

    u = zeros(nInputs,nSteps);
    % w = noise(1:nSteps,ii)';
    w = noise(1:nSteps,nStates*(ii-1)+1:nStates*(ii))';
    
    for k = 1:nSteps

        u(:,k) = K_opt*x(ii,:,k)';
        % x(ii,:,k+1) = (A*x(ii,:,k)' + B*u(:,k))'; %+ w(:,k))';
        x(ii,:,k+1) = (At*x(ii,:,k)' + Bt*u(:,k)+ w(:,k))';

    end

end

%% DLQR model based
Q = W_x;
R = W_u;
K_lqr = dlqr(At,Bt,Q,R);

x_lqr = zeros(nRealization,nStates,nSteps+1);

for ii = 1:nRealization

    % Simulation
    x_lqr(ii,:,1) = x0';

    u_lqr = zeros(nInputs,nSteps);
    % w = noise(1:nSteps,ii)';
    w = noise(1:nSteps,nStates*(ii-1)+1:nStates*(ii))';

    for k = 1:nSteps

        u_lqr(:,k) = -K_lqr*x_lqr(ii,:,k)';
        x_lqr(ii,:,k+1) = (At*x_lqr(ii,:,k)' + Bt*u_lqr(:,k)+ w(:,k))';

    end

end



%% Data only


Qo = sdpvar(nInputs,nStates);
Po = sdpvar(nStates,nStates);
L = sdpvar(nStates+nInputs, nStates+nInputs);


A_lm = [Po-eye(nStates), Qo',zeros(nStates);
        Qo, zeros(nInputs),-Qo;
        zeros(nStates),-Qo', Po];
B_lm = zeros(nStates+nInputs, nStates);
C_lm = -Po+eye(nStates);

mat_lm = -[C_lm, zeros(nStates), zeros(nStates,nInputs), zeros(nStates);
         zeros(nStates)',Po, Qo',zeros(nStates);
          zeros(nStates,nInputs)',Qo, zeros(nInputs),-Qo;
         zeros(nStates)',zeros(nStates),-Qo', -Po];


tau_dt = sdpvar(1,1);
miu_dt =sdpvar(1,1);
gamma = sdpvar(1,1);


constraintO_1 = - miu_dt * padded_mat_dt >= 0.0;
constraintO_2 = mat_lm - tau_dt*padded_mat_dt >=0;
constraintO_4 = tau_dt>=0.00000000001;
constraintO_5 = L - [0 , Qo;
                    Qo', Po]>=0;
constraintO_6 = trace(W_x*Po) +trace(W_u*L)<=gamma;
constraintO_7 = Po>=eye(nStates);

constraints_o_dt =   constraintO_1 + constraintO_2 + constraintO_4  + constraintO_5 + constraintO_6 + constraintO_7 ;
cost_o_dt = gamma;

ops = sdpsettings('solver','mosek','verbose',0);
diagnostics = optimize(constraints_o_dt, cost_o_dt, ops);

if(diagnostics.problem == 0)

    P_opt_opt = value(Po);
    Q_opt_opt = value(Qo);
    tau_dt_opt = value(tau_dt);
    % lambda_opt = value(lambda);

else

    fprintf('Error! Solution is infeasible!\n\n');

end

K_opt_dt = Q_opt_opt*inv(P_opt_opt);

x_opt_dt = zeros(nRealization,nStates,nSteps+1);

for ii = 1:nRealization


    % Simulation
    x_opt_dt(ii,:,1) = x0';

    u_opt_dt = zeros(nInputs,nSteps);
    % w = noise(1:nSteps,ii)';
    w = noise(1:nSteps,nStates*(ii-1)+1:nStates*(ii))';

    for k = 1:nSteps

        u_opt_dt(:,k) = K_opt_dt*x(ii,:,k)';
        x_opt_dt(ii,:,k+1) = (At*x_opt_dt(ii,:,k)' + Bt*u_opt_dt(:,k)+ w(:,k))';

    end

end

%% Uncertain only 

Qo_us = sdpvar(nInputs,nStates);
Po_us = sdpvar(nStates,nStates);
L_us = sdpvar(nStates+nInputs, nStates+nInputs);

A_lm = [Po_us-eye(nStates), Qo_us',zeros(nStates);
        Qo_us, zeros(nInputs),-Qo_us;
        zeros(nStates),-Qo_us', Po_us];
B_lm = zeros(nStates+nInputs, nStates);
C_lm = -Po_us+eye(nStates);

mat_lm = -[C_lm, zeros(nStates), zeros(nStates,nInputs), zeros(nStates);
         zeros(nStates)',Po_us, Qo_us',zeros(nStates);
          zeros(nStates,nInputs)',Qo_us, zeros(nInputs),-Qo_us;
         zeros(nStates)',zeros(nStates),-Qo_us', -Po_us];


tau_dt = sdpvar(1,1);
tau_us = sdpvar(1,1);
miu_us =sdpvar(1,1);
miu_dt =sdpvar(1,1);
gamma_us = sdpvar(1,1);


% constraintO_1_us = -miu_us * padded_mat_us >= 0.0;
constraintO_2_us = mat_lm - tau_us*padded_mat_us >=0;
constraintO_3_us = tau_us>=0.00000000001;
constraintO_5_us = L_us - [0 , Qo_us;
                    Qo_us', Po_us]>=0;
constraintO_6_us = trace(W_x*Po_us) +trace(W_u*L_us)<=gamma_us;
constraintO_7_us = Po_us>=eye(nStates);

constraints_o_us =   constraintO_2_us + constraintO_3_us  + constraintO_5_us + constraintO_6_us + constraintO_7_us ;
cost_o_us = gamma_us;

ops = sdpsettings('solver','mosek','verbose',0);
diagnostics = optimize(constraints_o_us, cost_o_us, ops);

if(diagnostics.problem == 0)

    P_opt_us = value(Po_us);
    Q_opt_us = value(Qo_us);
    % lambda_opt = value(lambda);

else

    fprintf('Error! Solution is infeasible!\n\n');

end

K_opt_uc = Q_opt_us*inv(P_opt_us);  

x_opt_uc = zeros(nRealization,nStates,nSteps+1);

for ii = 1:nRealization


    % Simulation
    x_opt_uc(ii,:,1) = x0';

    u_opt_uc = zeros(nInputs,nSteps);
    % w = noise(1:nSteps,ii)';
    w = noise(1:nSteps,nStates*(ii-1)+1:nStates*(ii))';

    for k = 1:nSteps

        u_opt_uc(:,k) = K_opt_uc*x(ii,:,k)';
        x_opt_uc(ii,:,k+1) = (At*x_opt_uc(ii,:,k)' + Bt*u_opt_uc(:,k)+ w(:,k))';
    end

end
%% plot 

cost_dt = zeros(nSteps,1);
cost_uc = zeros(nSteps,1);
cost_lqr = zeros(nSteps,1);
cost_uc_dt = zeros(nSteps,1);


for k=1:nSteps
    cost_temp_uc = 0;
    cost_temp_dt = 0;
    cost_temp_uc_dt = 0;
    cost_temp_lqr = 0;

    for j=1:k
        cost_temp_lqr = cost_temp_lqr + (x_lqr(ii,:,j)*W_x*x_lqr(ii,:,j)' + u_lqr(:,j)'*W_u*u_lqr(:,j));
        cost_temp_dt = cost_temp_dt + (x_opt_dt(ii,:,j)*W_x*x_opt_dt(ii,:,j)' + u_opt_dt(:,j)'*W_u*u_opt_dt(:,j));
        cost_temp_uc = cost_temp_uc + (x_opt_uc(ii,:,j)*W_x*x_opt_uc(ii,:,j)' + u_opt_uc(:,j)'*W_u*u_opt_uc(:,j));
        cost_temp_uc_dt = cost_temp_uc_dt + (x(ii,:,j)*W_x*x(ii,:,j)' + u(:,j)'*W_u*u(:,j));
    end
   cost_lqr(k) = cost_temp_lqr/k;
   cost_uc(k) = cost_temp_uc/k;
   cost_dt(k) = cost_temp_dt/k;
   cost_uc_dt(k) = cost_temp_uc_dt/k;
end



% Create the plot
figure;
h = plot(t, cost_uc, 'b-', t, cost_dt, 'r--', t, cost_uc_dt, 'g:', t, cost_lqr, 'm-.');
set(h, 'LineWidth', 3); % Set line width for all lines
xlabel('Steps');
ylabel('Average Cost');
% title('Costs vs. Steps');
legend('Uncertainty', 'Data', 'Physics-informed', 'LQR Model Based');
grid on;

% Adjust the appearance
set(gca, 'yscale', 'log'); % Use logarithmic scale for the y-axis

ax = gca;
ax.FontSize = 24; 

x_1_tf = reshape(x(ii,1,:),[nSteps+1,1]);
x_2_tf = reshape(x(ii,2,:),[nSteps+1,1]);
x_3_tf = reshape(x(ii,3,:),[nSteps+1,1]);

x_1_lqr = reshape(x_lqr(ii,1,:),[nSteps+1,1]);
x_2_lqr = reshape(x_lqr(ii,2,:),[nSteps+1,1]);
x_3_lqr = reshape(x_lqr(ii,3,:),[nSteps+1,1]);

x_1_dt = reshape(x_opt_dt(ii,1,:),[nSteps+1,1]);
x_2_dt = reshape(x_opt_dt(ii,2,:),[nSteps+1,1]);
x_3_dt = reshape(x_opt_dt(ii,3,:),[nSteps+1,1]);

x_1_uc = reshape(x_opt_uc(ii,1,:),[nSteps+1,1]);
x_2_uc = reshape(x_opt_uc(ii,2,:),[nSteps+1,1]);
x_3_uc = reshape(x_opt_uc(ii,3,:),[nSteps+1,1]);

% Create a new figure and set its position
figure;
set(gcf, 'Position', [100, 100, 800, 800]);

% Create the subplots
subplot(3, 1, 1);

h = plot(t_, x_1_uc, 'b-', t_, x_1_dt, 'r--', t_, x_1_tf, 'g:', t_, x_1_lqr, 'm-.');
set(h, 'LineWidth', 3); % Set line width for all lines
xlabel('Steps');
ylabel('V1');
legend('Uncertainty', 'Data', 'Physics-informed', 'LQR Model Based');
grid on;

ax = gca;
ax.FontSize = 24; 

subplot(3, 1, 2);

h = plot(t_, x_2_uc, 'b-', t_, x_2_dt, 'r--', t_, x_2_tf, 'g:', t_, x_2_lqr, 'm-.');
set(h, 'LineWidth', 3); % Set line width for all lines
xlabel('Steps');
ylabel('V2');
legend('Uncertainty', 'Data', 'Physics-informed', 'LQR Model Based');
grid on;

ax = gca;
ax.FontSize = 24; 


subplot(3, 1, 3);

h = plot(t_, x_3_uc, 'b-', t_, x_3_dt, 'r--', t_, x_3_tf, 'g:', t_, x_3_lqr, 'm-.');
set(h, 'LineWidth', 3); % Set line width for all lines
xlabel('Steps');
ylabel('SoC');
legend('Uncertainty', 'Data', 'Physics-informed', 'LQR Model Based');
grid on;

ax = gca;
ax.FontSize = 24; 

% Create the plot
figure;
h = plot(t, u_opt_uc, 'b-', t, u_opt_dt, 'r--', t, u, 'g:', t, u_lqr, 'm-.');
set(h, 'LineWidth', 3); % Set line width for all lines
xlabel('Steps');
ylabel('Control input');
% title('Costs vs. Steps');
legend('Uncertainty', 'Data', 'Physics-informed', 'LQR Model Based');
grid on;

ax = gca;
ax.FontSize = 24; 








