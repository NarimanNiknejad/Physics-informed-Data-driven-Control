clc;
clear;
close all;

%% Define general parameters
nStates = 2;
nInputs = 1;
nRealization = 100;
nSample = 100000;

lambda = 0.9;

%% Time sequence
T0 = 0;
Ts = 0.01;
Tf = 20;

t = T0:Ts:Tf;

nSteps = numel(t);
%% Noise generation
mu = [0 0];
varNoise = 0.0005;
coVarNoise = 0;
coVarMatrix = [  varNoise   coVarNoise
                coVarNoise   varNoise  ];
sigma = coVarMatrix;

noise = noise_generator(mu,sigma,nStates,nSample,nRealization);
%% Drone Simulation Model 
m = 1.380;
g = 9.81;

% Z - axis
Az = [0.0, 1.0;
      0.0, 0.0];
Bz = [0.0;
     1 / (m)];
Cz = [1, 0];

% Finding the discrete-time dynamics of the drone 
time_step = 0.01;
sys = ss(Az,Bz,Cz,[]);
dis = c2d(sys,time_step,'zoh');
% Az_d = dis.A;
% Bz_d = dis.B;
% Cz_d = dis.C;

At = [0.95, 0.1;
    0, 0.95];
Bt = [0.0000;
    0.0072];

%% Define the parameters related to the actual system and the nominal system

% The difference 
ep_AB = 0.001;
delta_A =  ep_AB*eye(nStates);
delta_B = ep_AB*[0;1];

As = At - delta_A;
Bs = Bt - delta_B;


%% Define parameters for the physics-informed set

SourceSys = [As, Bs]';
A_us = eye(nStates+nInputs);
B_us = -[As, Bs]';
C_us = SourceSys'*SourceSys - ep_AB^2*eye(nStates);

mat_us = -[C_us, B_us';
           B_us, A_us ];

padded_mat_us = padarray(mat_us, [2 2], 0, 'post');
%% Define parameters for collected data

% First we need to collect some data from the actual system
N = nStates + 1;

X = zeros(nStates,N+1);
X(:,1) = [0 1]';
U0 = 5*randn(nInputs,N);
W0 = noise(1:N,1:2)';

for kk = 1:N

    X(:,kk+1) = At*X(:,kk) + Bt*U0(:,kk) + W0(:,kk);
    % X(:,kk+1) = At*X(:,kk) + Bt*U0(:,kk); %no noise

end

X0 = X(:,1:N);
X1 = X(:,2:N+1);

data = [X0; U0];
data_use = data;
% Check rank
if(rank(data)<min(size(data)))

    errordlg('Data is not full row rank!');

end

% Now we define the parameters for the data-driven uncertainty
A_dt = data_use*data_use';
B_dt = -data_use*X1';
C_dt = -sigma + X1*X1';

mat_dt = -[  C_dt, B_dt';
             B_dt, A_dt ];

padded_mat_dt = padarray(mat_dt, [2 2], 0, 'post');

%% Define parameters for Ph-DD safe algorithm

Qs = sdpvar(nInputs,nStates);
Ps = sdpvar(nStates,nStates);
% lambda = sdpvar(1,1);
P_safe =  [0.7, 0.010;
0.010, 0.10];


A_lm = [Ps, Qs',zeros(nStates);
        Qs, zeros(nInputs),-Qs;
        zeros(nStates),-Qs', Ps];
B_lm = zeros(nStates+nInputs, nStates);
C_lm = -lambda^2*Ps;

mat_lm = -[C_lm, zeros(nStates), zeros(nStates,nInputs), zeros(nStates);
         zeros(nStates)',Ps, Qs',zeros(nStates);
          zeros(nStates,nInputs)',Qs, zeros(nInputs),-Qs;
         zeros(nStates)',zeros(nStates),-Qs', -Ps];


tau_dt = sdpvar(1,1);
tau_us = sdpvar(1,1);
miu_us =sdpvar(1,1);
miu_dt =sdpvar(1,1);


constraintS_1 = -miu_us * padded_mat_us - miu_dt * padded_mat_dt >= 0.0;
constraintS_2 = mat_lm - tau_us*padded_mat_us - tau_dt*padded_mat_dt >=0;
constraintS_3 = tau_us>=0.000000001;
constraintS_4 = tau_dt>=0.000000001;
constraintS_5 = Ps>=P_safe;


constraints_s =   constraintS_1 + constraintS_2 + constraintS_3 + constraintS_4  + constraintS_5;
cost_s = trace(Ps);

ops = sdpsettings('solver','mosek','verbose',0);
diagnostics_s = optimize(constraints_s, cost_s, ops);

if(diagnostics_s.problem == 0)

    P_safe_opt = value(Ps);
    Q_safe_opt = value(Qs);
    tau_us_opt = value(tau_us);
    tau_dt_opt = value(tau_dt);
    % lambda_opt = value(lambda);

else

    fprintf('Error! Solution is infeasible!\n\n');

end

K_safe = Q_safe_opt*inv(P_safe_opt);  



%% Main loop
x = zeros(nRealization,nStates,nSteps+1);
x0 = [-0.84, 2.32]';

for ii = 1:nRealization


    % Simulation
    x(ii,:,1) = x0';

    u = zeros(nInputs,nSteps);
    % w = noise(1:nSteps,ii)';
    w = noise(1:nSteps,nStates*(ii-1)+1:nStates*(ii))';
    

    for k = 1:nSteps

        x_now = x(ii,:,k);

        alpha_opt = 1;
        u(:,k) = (K_safe)*x(ii,:,k)';
        % x(ii,:,k+1) = (A*x(ii,:,k)' + B*u(:,k))'; %+ w(:,k))';
        x(ii,:,k+1) = (At*x(ii,:,k)' + Bt*u(:,k)+ w(:,k))';

    end

end

%% Plot results
% Phase plane
colors = jet(nRealization);

fig1 = figure(1);
% subplot(1,2,1);
ops.plot.edgeColor = 'g';
ops.plot.wirestyle = '--';
xx = sdpvar(nStates,1);


alpha(0.0);

for i = 1:30

    plot(xx'*P_safe_opt*xx<=lambda^(i-1),[],'b');
    hold on;

end

% 
plot(xx'*P_safe*xx<=1,[],'r');
alpha(0.1);


plot(x0(1),x0(2),'sm','MarkerSize',15,'MarkerFaceColor','m');
for ii = 1:nRealization

    x_ = reshape(x(ii,1,:),[nSteps+1,1]);
    y_ = reshape(x(ii,2,:),[nSteps+1,1]);
    plot(x_,y_,'Color', colors(ii, :),'LineWidth',1,'MarkerSize',12);

end
grid on;
axis equal;
xlabel('x_1 (t)','FontSize',24);
ylabel('x_2 (t)','FontSize',24);
legend('','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','', ...
       'Initial point','FontSize',24,'TextColor','k');
% legend('boxoff');
ax = gca;
ax.FontSize = 24; 

