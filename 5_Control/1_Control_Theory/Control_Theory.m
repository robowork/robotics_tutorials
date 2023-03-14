%% Helper to define transfer functions more quickly
s = tf([1 0],1)

%% First-Order system
omega_0 = 1

firstorder = 1/((s/omega_0) + 1)

figure;
pzmap(firstorder);

figure;
step(firstorder);

figure;
bode(firstorder);

%% Second-Order system
omega_0 = 1;
zeta = 0.5;

secondorder = 1/((s/omega_0)^2 + 2*zeta*(s/omega_0)+1)

figure;
pzmap(secondorder);

figure;
step(secondorder);

figure;
bode(secondorder);

%% Integrator
integrator = 1/s

figure;
bode(integrator);

%% Differentiator
differentiator = 1/integrator

figure;
bode(differentiator);

%% Closed-Loop (Unit Gain) Stability

figure;
hold on;
nyquist(firstorder);
nyquist(secondorder);
legend('firstorder','secondorder');
hold off;


unstablesys = (0.04*s^2+0.12*s)/(s^3+0.1*s^2+0.8*s-0.4)

format long;  %to show more decimals for poles
pole(unstablesys)
format short;  %back to regular formatting

figure;
pzmap(unstablesys);
legend('unstablesys');

figure;
step(unstablesys);
legend('unstablesys');

% Note: open-loop system is unstable (1 unstable pole), regular Bode margins approach cannot be used
figure;
bode(unstablesys);
legend('unstablesys');

% Note: 0 encirclements of -1, therefore closed-loop system will have the same number of unstable poles as the open-loop one (we know it has 1)
figure;
nyquist(unstablesys);
legend('unstablesys');


f16 = (3.553e-15*s^4-0.1642*s^3-0.1243*s^2-0.00161*s+9.121e-17)/(s^5+1.825*s^4+2.941*s^3+0.03508*s^2+0.01522*s-1.245e-15)

format long;  %to show more decimals for poles
pole(f16)
format short;  %back to regular formatting

figure;
pzmap(f16);
legend('f16');

% Note: open-loop system is unstable (1 unstable pole), regular Bode margins approach cannot be used
figure;
bode(f16);
legend('f16');

% Note: 2 clockwise encirclements of -1, therefore closed-loop system will have 2 more unstable poles than the open-loop one (we know it has 1)
figure;
nyquist(f16);
legend('f16');

%% Minimum Phase / Non-Minimum Phase systems
G_minphase = (s+2)/(s^2+3*s+1);
G_inputdelay = tf(1,1,'InputDelay',1) * G_minphase;
G_nonminphase = (s-2)/(s^2+3*s+1);

opts = bodeoptions('cstprefs');
opts.PhaseWrapping = 'on';

hold on;
bode(G_minphase,opts);
bode(G_inputdelay,opts);
bode(G_nonminphase,opts);
legend('G\_minphase','G\_inputdelay','G\_nonminphase');
hold off;

%% PID Control
I_yy = 0.1;

pitch = 1/(I_yy*s^2)

K_p = 1;
K_d = 1;
K_i = 0.1;

controller_P = K_p
controller_PD = K_p + K_d*s
controller_PI = K_p + K_i/s
controller_PID = K_p + K_d*s + K_i/s

pitch_P = feedback(controller_P * pitch, 1)
pitch_PD = feedback(controller_PD * pitch, 1)
pitch_PI = feedback(controller_PI * pitch, 1)
pitch_PID = feedback(controller_PID * pitch, 1)

hold on;
step(pitch,30);
step(pitch_P,30);
step(pitch_PD,30);
step(pitch_PI,30);
step(pitch_PID,30);
set(gca,'ylim',[-10 10]);
legend('open','P','PD','PI','PID');
hold off;

%% State-Space
m = 1;
k = 1;
b = 1;

A = [0   1;
    -k/m -b/m]
B = [0; 
     1/m]
C = [1 0;
     0 1]
D = [0;
     0]

harmonic_oscillator = ss(A,B,C,D)

figure;
step(harmonic_oscillator)
legend('harmonic\_oscillator');

figure;
impulse(harmonic_oscillator)
legend('harmonic\_oscillator');

[numerator,denominator] = ss2tf(A,B,C,D)
output1_tf = tf(numerator(1,:),denominator)
output2_tf = tf(numerator(2,:),denominator)

% Eigenvectors,Eigenvalues
[eigenvectors,eigenvalues] = eig(harmonic_oscillator.A)
eigenvector1 = eigenvectors(:,1)
eigenvector2 = eigenvectors(:,2)
eigenvalue1 = eigenvalues(1,1)
eigenvalue2 = eigenvalues(2,2)
harmonic_oscillator.A*[eigenvector1 eigenvector2] - [eigenvector1 eigenvector2]*diag([eigenvalue1 eigenvalue2])  %check equality (should be close to 0 matrix)
pole(output1_tf) - [eigenvalue1;eigenvalue2]  %check equality (should be close to 0 column)
pole(output2_tf) - [eigenvalue1;eigenvalue2]  %check equality (should be close to 0 column)

% Stability
real([eigenvalue1 eigenvalue2])  %check stability (should all be <0)

[harmonic_oscillator_rows,harmonic_oscillator_cols] = size(harmonic_oscillator.A);

% Controlability
C_ctrb = ctrb(harmonic_oscillator)
rank(C_ctrb) == harmonic_oscillator_cols  % controllable

% Note: How to check for stabilizability IF system is not controllable
% [harmonic_oscillator_stable,harmonic_oscillator_unstable] = stabsep(harmonic_oscillator)
% [dynamics_ss_unstable_rows,dynamics_ss_unstable_cols] = size(harmonic_oscillator.A);
% C_stab = ctrb(harmonic_oscillator_unstable)
% rank(C_stab) == dynamics_ss_unstable_cols  % stabilizable ?

% Observability
O_obsv = obsv(harmonic_oscillator)
rank(O_obsv) == harmonic_oscillator_cols  % observable
