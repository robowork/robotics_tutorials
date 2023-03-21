%% Helper to define transfer functions more quickly
s = tf([1 0],1)

%% System definition

% Pitch dynamics
I_yy = 1e-4;
pitch_tf = 1/(I_yy*s^2)
pitch_ss = ssform(idss(pitch_tf),'Form','canonical');

% Use PIDTool to find working PD gains
%pidtool(pitch_ss)  

K_p = 0.05;
K_d = 0.0025;
pitch_controller_PD = K_p + K_d*s
pitch_closedloop_ss = minreal( feedback(pitch_controller_PD * pitch_ss, 1) )

figure;
step(pitch_closedloop_ss)
legend('pitch dynamics (PD-controlled closedloop) step');

% Longitudinal dynamics
dx_tf = -9.81/(s^2)
dx_ss = ssform(idss(dx_tf),'Form','canonical')

dynamics_ss = series(pitch_closedloop_ss,dx_ss)

figure;
step(dynamics_ss)
legend('longitudinal dynamics step');

% Tests
[dynamics_ss_rows,dynamics_ss_cols] = size(dynamics_ss.A);

C_ctrb = ctrb(dynamics_ss);
rank(C_ctrb) == dynamics_ss_cols  % controllable

O_obsv = obsv(dynamics_ss);
rank(O_obsv) == dynamics_ss_cols  % observable

%% LQR - Continuous -vs- Continuous with Input & Input Rate Constraints

Q = diag([10 0.1 0.0001 0])
R = diag([1])

num_states = 4;
num_outputs = 1;
num_inputs = 1;

K_lqr = lqr(dynamics_ss,Q,R)

dynamics_ss_lqr_Aclosedloop = dynamics_ss.A - dynamics_ss.B*K_lqr;
dynamics_ss_lqr_closedloop = ss(dynamics_ss_lqr_Aclosedloop,dynamics_ss.B,dynamics_ss.C,dynamics_ss.D)

dcgain_compensator = tf([1/dcgain(dynamics_ss_lqr_closedloop)],[1]);
K_r = cell2mat(dcgain_compensator.Numerator);

dynamics_ss_lqr_closedloop = series(K_r, dynamics_ss_lqr_closedloop)


dt = 0.01;
t = linspace(0,15,15/dt);

figure;
step(dynamics_ss_lqr_closedloop,t)
legend('LQR continuous');


system_sim = dynamics_ss;

x_0 = zeros(num_states,1);
x = x_0;

x_t = [];
y_t = [];
u_t = [];
for i = 1:length(t)
    t_i = t(i);
    x_t = [x_t; x'];
    
    r = 1;
    u = K_r * r - K_lqr * x;
    
    
    % Ad-hoc enforcement of constraints
    u = max(min(u,deg2rad(10)),deg2rad(-10)); % constraint
    if isempty(u_t)
        u_prev = zeros(num_inputs,1);
    else
        u_prev = u_t(end,1)';
    end
    dudt = (u-u_prev)/dt;
    dudt = max(min(dudt,deg2rad(10)),deg2rad(-10)); % rate constraint
    u = u_prev + dudt * dt;
    
    
    x_dot = system_sim.A * x + system_sim.B * u;
    y = system_sim.C * x + system_sim.D * u;

    x = x + x_dot * dt;
    
    y_t = [y_t; y'];
    u_t = [u_t; u'];
end

figure;
yyaxis left;
plot(t,y_t(:,1));
yyaxis right;
plot(t,u_t(:,1));
legend('LQR continuous(& adhoc constraints) output','LQR continuous (& adhoc constraints) input');

%% LQR - Discrete -vs- Discrete with Input & Input Rate Constraints

Q = diag([10 0.1 0.0001 0])
R = diag([1])
N = zeros(4,1)

num_states = 4;
num_outputs = 1;
num_inputs = 1;

Ts = 0.1;
dynamics_ss_discrete = c2d(dynamics_ss,Ts);

K_dlqr = dlqr(dynamics_ss_discrete.A,dynamics_ss_discrete.B,Q,R,N)

dynamics_ss_discrete_lqr_Aclosedloop = dynamics_ss_discrete.A - dynamics_ss_discrete.B*K_dlqr;
dynamics_ss_discrete_lqr_closedloop = ss(dynamics_ss_discrete_lqr_Aclosedloop,dynamics_ss_discrete.B,dynamics_ss_discrete.C,dynamics_ss_discrete.D,Ts)

dcgain_compensator_discrete = tf([1/dcgain(dynamics_ss_discrete_lqr_closedloop)],[1]);
K_dr = cell2mat(dcgain_compensator_discrete.Numerator);

dynamics_ss_discrete_lqr_closedloop = series(K_dr, dynamics_ss_discrete_lqr_closedloop)


Ts = Ts;
t = [0:Ts:15-Ts];

figure;
step(dynamics_ss_discrete_lqr_closedloop,t)
legend('LQR discrete (unconstrained)');


system_sim = dynamics_ss_discrete;

x_0 = zeros(num_states,1);
x = x_0;

x_t = [];
y_t = [];
u_t = [];
for i = 1:length(t)
    t_i = t(i);
    x_t = [x_t; x'];
    
    r = 1;
    u = K_dr * r - K_dlqr * x;
    
    
    % Ad-hoc enforcement of constraints
    u = max(min(u,deg2rad(10)),deg2rad(-10)); % constraint
    if isempty(u_t)
        u_prev = zeros(num_inputs,1); 
    else
        u_prev = u_t(end,1)';
    end
    dudt = (u-u_prev)/Ts;
    dudt = max(min(dudt,deg2rad(10)),deg2rad(-10)); % rate constraint
    u = u_prev + dudt * Ts;
    
    
    x_new = system_sim.A * x + system_sim.B * u;
    y = system_sim.C * x + system_sim.D * u;

    x = x_new;
    
    y_t = [y_t; y'];
    u_t = [u_t; u'];
end

figure;
yyaxis left;
stairs(t,y_t(:,1));
yyaxis right;
stairs(t,u_t(:,1));
legend('LQR discrete (& adhoc constraints) output','LQR discrete (& adhoc constraints) input');

%% MPC (Receding Horizon Unconstrained Optimization)
clc;

Q = diag([1]);  % Outputs Q
R = diag([1]);  % Inputs R

N_horizon = 10;

num_states = 4;
num_outputs = 1;
num_inputs = 1;

Y_horizon.W = zeros(N_horizon*num_outputs , num_states);
Y_horizon.Z = zeros(N_horizon*num_outputs , N_horizon);
for j = 1:N_horizon
    Y_horizon.W(j:j-1+num_outputs,:) = dynamics_ss_discrete.C * (dynamics_ss_discrete.A^j);
    
    for k = 1:j
        Y_horizon.Z(j:j-1+num_outputs,k) = dynamics_ss_discrete.C * (dynamics_ss_discrete.A^(j-k)) * dynamics_ss_discrete.B;
    end
end

r = ones(N_horizon*num_outputs, 1);

x_0 = 0.5 * ones(num_states,1);

DU = 0.1 * ones(N_horizon*num_outputs, 1);

Q_blkdiag = [];
R_blkdiag = [];
for j = 1:N_horizon
    Q_blkdiag = blkdiag(Q_blkdiag,Q);
    R_blkdiag = blkdiag(R_blkdiag,R);
end

format long;

% Dimensionality sanity check
J = 0.5 * (r - (Y_horizon.W * x_0 + Y_horizon.Z * DU))' * Q_blkdiag * (r - (Y_horizon.W * x_0 + Y_horizon.Z * DU)) + ...
    0.5 * (DU)' * R_blkdiag * DU

% Quadratic breakdown ( z'*Q*z for z=a+bu => (bu)'*Q*(bu)+2*a'*(bu) => u'*(b'Q*)*u+2*(a'b)*u) sanity check
J = 0.5 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * (r - (Y_horizon.W * x_0)) + ...
    0.5 * (-(Y_horizon.Z * DU))' * Q_blkdiag * (-(Y_horizon.Z * DU)) + ...
    0.5 * 2 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * (-(Y_horizon.Z * DU)) + ...
    0.5 * (DU)' * R_blkdiag * DU

% Quadratic refactoring into: ct + (1/2)*DU'*A*DU - b'*DU sanity check
J = (0.5 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * (r - (Y_horizon.W * x_0))) + ...  % constant term
    0.5 * DU' * ((-Y_horizon.Z)'*Q_blkdiag*(-Y_horizon.Z) + R_blkdiag) * DU + ...  % DU quadratic term
    ((0.5 * 2 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * -(Y_horizon.Z))')' * DU  % DU linear term

% Final sanity check
ct_quad = (0.5 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * (r - (Y_horizon.W * x_0)));
A_quad = (-Y_horizon.Z)'*Q_blkdiag*(-Y_horizon.Z) + R_blkdiag;
b_quad = (0.5 * 2 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * -(Y_horizon.Z))';
J = ct_quad + ...
    (1/2) * DU' * A_quad * DU + ...
    b_quad' * DU

format short;

% Sample solution of unconstrained quadratic problem for 1 prediction horizon (https://sites.math.washington.edu/~burke/crs/408/notes/Math408_W2020/quad_opt.pdf)
DU_opt = -inv(A_quad) * b_quad


% Simulate
system_sim = dynamics_ss_discrete;

x_0 = zeros(num_states,1);
x = x_0;

x_t = [];
y_t = [];
u_t = [];
for i = 1:length(t)
    t_i = t(i);
    x_t = [x_t; x'];
    
    x_0 = x;
    A_quad = (-Y_horizon.Z)'*Q_blkdiag*(-Y_horizon.Z) + R_blkdiag;
    b_quad = (0.5 * 2 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * -(Y_horizon.Z))';
    DU_opt = -inv(A_quad) * b_quad;

    if isempty(u_t)
        u_prev = zeros(num_inputs,1);
    else
        u_prev = u_t(end,1)';
    end
    u = DU_opt(1);  % keep only first control move in optimal sequence
        
    
    % Ad-hoc enforcement of constraints
    u = max(min(u,deg2rad(10)),deg2rad(-10)); % constraint
    if isempty(u_t)
        u_prev = zeros(num_inputs,1); 
    else
        u_prev = u_t(end,1)';
    end
    dudt = (u-u_prev)/Ts;
    dudt = max(min(dudt,deg2rad(10)),deg2rad(-10)); % rate constraint
    u = u_prev + dudt * Ts;
    
    
    x_new = system_sim.A * x + system_sim.B * u;
    y = system_sim.C * x + system_sim.D * u;

    x = x_new;
    
    y_t = [y_t; y'];
    u_t = [u_t; u'];
end

figure;
yyaxis left;
stairs(t,y_t(:,1));
yyaxis right;
stairs(t,u_t(:,1));
legend('MPC (unconstrained optimization & adhoc constraints) output','MPC (unconstrained optimization & adhoc constraints) input');

%% MPC with Constrained Optimization
clc;

% LMI constraints

DU_lb = deg2rad(-10);
DU_ub = deg2rad(10);
DU_lb_stacked = [];
DU_ub_stacked = [];

DU_lhs = [-1; ...
          1];
DU_rhs = [deg2rad(10); ...
          deg2rad(10)];
DU_lhs_stacked = [];
DU_rhs_stacked = [];

U_rate_lhs = [-1; ...
              1];
U_rate_rhs = [deg2rad(10)*Ts; ...
              deg2rad(10)*Ts];
U_rate_lhs_stacked = zeros(N_horizon*2*num_inputs,N_horizon);
U_rate_rhs_stacked = [];

for j = 1:N_horizon
    DU_lb_stacked = [DU_lb_stacked;DU_lb];
    DU_ub_stacked = [DU_ub_stacked;DU_ub];
    
    DU_lhs_stacked = blkdiag(DU_lhs_stacked,DU_lhs);
    DU_rhs_stacked = [DU_rhs_stacked;DU_rhs];

%     for k = 1:j
%         for l = 1:num_outputs
%             U_rate_lhs_stacked(2*(j-1+l)-1,k) = -1;
%             U_rate_lhs_stacked(2*(j-1+l),k) = 1;
%         end
%     end
    for l = 1:num_outputs
        U_rate_lhs_stacked(2*(j-1+l)-1,j) = -1;
        U_rate_lhs_stacked(2*(j-1+l),j) = 1;
    end
    k = j-1;
    if ~(k <= 0)
        U_rate_lhs_stacked(2*(j-1+l)-1,k) = -(-1);
        U_rate_lhs_stacked(2*(j-1+l),k) = -(1);
    end
    U_rate_rhs_stacked = [U_rate_rhs_stacked;U_rate_rhs];
end

% Simulate
system_sim = dynamics_ss_discrete;

x_0 = zeros(num_states,1);
x = x_0;

x_t = [];
y_t = [];
u_t = [];
for i = 1:length(t)
    t_i = t(i);
    x_t = [x_t; x'];
    
    x_0 = x;
    A_quad = (-Y_horizon.Z)'*Q_blkdiag*(-Y_horizon.Z) + R_blkdiag;
    b_quad = (0.5 * 2 * (r - (Y_horizon.W * x_0))' * Q_blkdiag * -(Y_horizon.Z))';   
    % DU_opt = quadprog(A_quad,b_quad)  % unconstrained
    % DU_opt = quadprog(A_quad,b_quad,[],[],[],[],DU_lb_stacked,DU_ub_stacked)  % alternative (simpler) way to specify lower/upper bounds in Matlab's quadprog
    % DU_opt = quadprog(A_quad,b_quad,DU_lhs_stacked,DU_rhs_stacked)  % only input magnitude constraints
    % DU_opt = quadprog(A_quad,b_quad,U_rate_lhs_stacked,U_rate_rhs_stacked)  % only input rate constraints
    DU_opt = quadprog(A_quad,b_quad,[DU_lhs_stacked;U_rate_lhs_stacked],[DU_rhs_stacked;U_rate_rhs_stacked])  % input magnitude and input rate constraints

    if isempty(u_t)
        u_prev = zeros(num_inputs,1);
    else
        u_prev = u_t(end,1)';
    end
    u = DU_opt(1);  % keep only first control move in optimal sequence
        
    
    % Ad-hoc enforcement of constraints
    u = max(min(u,deg2rad(10)),deg2rad(-10)); % constraint
    if isempty(u_t)
        u_prev = zeros(num_inputs,1); 
    else
        u_prev = u_t(end,1)';
    end
    dudt = (u-u_prev)/Ts;
    dudt = max(min(dudt,deg2rad(10)),deg2rad(-10)); % rate constraint
    u = u_prev + dudt * Ts;
    
    
    x_new = system_sim.A * x + system_sim.B * u;
    y = system_sim.C * x + system_sim.D * u;

    x = x_new;
    
    y_t = [y_t; y'];
    u_t = [u_t; u'];
end

figure;
yyaxis left;
stairs(t,y_t(:,1));
yyaxis right;
stairs(t,u_t(:,1));
legend('MPC (constrained optimization & adhoc constraints) output','MPC (constrained optimization & adhoc constraints) input');