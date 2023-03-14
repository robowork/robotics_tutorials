%% Helper to define transfer functions more quickly
s = tf([1 0],1)

%% System definition

% Actuator dynamics
communication_transportdelay = 0.1;

actuator_omega_0 = 10/2*pi
actuator_tf = 1/((s/actuator_omega_0) + 1)
actuator_ss = ssform(idss(actuator_tf,'InputDelay',communication_transportdelay),'Form','canonical')
figure;
step(actuator_ss)
legend('actuator step');

% Pitch dynamics
I_yy = 1e-4;
pitch_tf = 1/(I_yy*s^2)
pitch_ss = ssform(idss(pitch_tf),'Form','canonical');
pitch_with_actuator_ss = series(pitch_ss,actuator_ss)

% Use PIDTool to find working PD gains
%pidtool(pitch_with_actuator_ss)  

K_p = 5e-5;
K_d = 10*K_p;
pitch_controller_PD = K_p + K_d*s
pitch_closedloop_ss = minreal( feedback(pitch_controller_PD * pitch_with_actuator_ss, 1) )

figure;
step(pitch_closedloop_ss)
legend('pitch dynamics (PD-controlled closedloop) step');

% Longitudinal dynamics
dx_tf = -9.81/(s^2)
dx_ss = ssform(idss(dx_tf),'Form','canonical')

dynamics_ss = pade( series(pitch_closedloop_ss,dx_ss) )

figure;
step(dynamics_ss)
legend('longitudinal dynamics (Pade-approximated) step');

% Tests
[dynamics_ss_rows,dynamics_ss_cols] = size(dynamics_ss.A);

C_ctrb = ctrb(dynamics_ss);
rank(C_ctrb) == dynamics_ss_cols  % controllable

O_obsv = obsv(dynamics_ss);
rank(O_obsv) == dynamics_ss_cols  % observable

%% Pole Placement

pole(dynamics_ss)

pole_new = [-0.01 -0.05 -0.1+0.1i -0.1-0.1i -0.1+0.2i -0.1-0.2i] 

K_place = place(dynamics_ss.A,dynamics_ss.B,pole_new)

dynamics_ss_place_Aclosedloop = dynamics_ss.A - dynamics_ss.B*K_place;
dynamics_ss_place_closedloop = ss(dynamics_ss_place_Aclosedloop,dynamics_ss.B,dynamics_ss.C,dynamics_ss.D)

dcgain_compensator = tf([1/dcgain(dynamics_ss_place_closedloop)],[1]);

dynamics_ss_place_closedloop = series(dcgain_compensator , dynamics_ss_place_closedloop)

figure;
step(dynamics_ss_place_closedloop)
legend('PolePlacement step');

%% LQR

Q = diag([1 0.1 0 0 0 0])
R = diag([1])

K_lqr = lqr(dynamics_ss,Q,R)

dynamics_ss_lqr_Aclosedloop = dynamics_ss.A - dynamics_ss.B*K_lqr;
dynamics_ss_lqr_closedloop = ss(dynamics_ss_lqr_Aclosedloop,dynamics_ss.B,dynamics_ss.C,dynamics_ss.D)

dcgain_compensator = tf([1/dcgain(dynamics_ss_lqr_closedloop)],[1]);

dynamics_ss_lqr_closedloop = series(dcgain_compensator , dynamics_ss_lqr_closedloop)

figure;
step(dynamics_ss_lqr_closedloop)
legend('LQR');

%% Extra: Helper system to recover input of dx_ss (pitch angle) block in series system (pitch_closedloop_ss->dx_ss) output for inspection

feedthrough_ss = ssform(idss(ss(0,0,0,1)),'Form','canonical')
dx_with_pitchfeedthrough = append(dx_ss,feedthrough_ss)
dynamics_with_pitchfeedthrough_ss = series(pitch_closedloop_ss,dx_with_pitchfeedthrough)

% Augment LQR controlled system to ouput a feedthrough of the dx_ss input (pitch angle) for inspection
dynamics_ss_lqr_Aclosedloop_augmented = [dynamics_ss_lqr_closedloop.A zeros(length(dynamics_ss_lqr_closedloop.A),length(dynamics_with_pitchfeedthrough_ss.A)-length(dynamics_ss_lqr_closedloop.A));
                                         dynamics_with_pitchfeedthrough_ss.A(length(dynamics_ss_lqr_closedloop.A)+1:end,:)]
dynamics_ss_lqr_Bclosedloop_augmented = [dynamics_ss_lqr_closedloop.B;
                                         dynamics_with_pitchfeedthrough_ss.B(length(dynamics_ss_lqr_closedloop.A)+1:end,1)]
dynamics_ss_lqr_Cclosedloop_augmented = [dynamics_ss_lqr_closedloop.C zeros(1,length(dynamics_with_pitchfeedthrough_ss.A)-length(dynamics_ss_lqr_closedloop.A));
                                         dynamics_with_pitchfeedthrough_ss.C(2,:)]
dynamics_ss_lqr_Dclosedloop_augmented = [dynamics_ss_lqr_closedloop.D];
                                     
dynamics_ss_lqr_closedloop_augmented = ss(dynamics_ss_lqr_Aclosedloop_augmented,dynamics_ss_lqr_Bclosedloop_augmented,dynamics_ss_lqr_Cclosedloop_augmented,dynamics_ss_lqr_Dclosedloop_augmented)

figure;
step(dynamics_ss_lqr_closedloop_augmented)
legend('LQR augmented system');