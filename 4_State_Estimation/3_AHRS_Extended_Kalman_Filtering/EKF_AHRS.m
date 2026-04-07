%% Generation of IMU & Magnetometer measurement data %%
clear all; clc;
g = 9.81;  %m/s^2
mag_G = 0.5;  %Gauss
sigma_pqr = deg2rad(30);  %+-15 deg/s
sigma_xyz = 0.1*g;  %+-10% of 9.81 m/s^2
sigma_mag = 0.1*mag_G;  %+-10% of mag_G

gyro_randwalk_enable = 0;

Ts = 0.01;
t = [0:Ts:10];

% Create some rotational time sequence (attention! these values cannot be used "as is", since pitch could go over vertical,
% but they can be used to generate equivalent SO(3) rotations and work with these onwards...) 
scalefactor = 0.5;
imu_yaw = wrapToPi( scalefactor * deg2rad(270)*sin(1*t) );
imu_pitch = wrapToPi( scalefactor * deg2rad(120)*sin(0.5*t) );
imu_roll = wrapToPi( scalefactor * deg2rad(90)*sin(2*t) );

figure;
plot(t,imu_yaw,t,imu_pitch,t,imu_roll);
legend('imu yaw','imu pitch','imu roll');

imu_R_SB = {};
for i=1:length(t)
  imu_R_SB{i} = eul2rotm([imu_yaw(i) imu_pitch(i) imu_roll(i)],'ZYX');
end

% Calculate (body) rotational rates and generate noise-corrupted measurements
imu_omega_B = {};
for i=1:length(t)
  if i == length(t)
      imu_omega_B{i} = imu_omega_B{i-1};
      break;
  end
  imu_R_SB_dot = (imu_R_SB{i+1}-imu_R_SB{i})/Ts;
  imu_omega_B{i} = imu_R_SB{i}' * imu_R_SB_dot;
  imu_omega_B{i}(1,1) = 0; imu_omega_B{i}(2,2) = 0; imu_omega_B{i}(3,3) = 0;  %force skew-symmetric 
end

imu_r = NaN*ones(size(t));
imu_q = NaN*ones(size(t));
imu_p = NaN*ones(size(t));
meas_r = NaN*ones(size(t));
meas_q = NaN*ones(size(t));
meas_p = NaN*ones(size(t));
randwalk_r = 0;
randwalk_q = 0;
randwalk_p = 0;
for i=1:length(t)
    imu_r(i) = -imu_omega_B{i}(1,2);
    imu_q(i) = imu_omega_B{i}(1,3);
    imu_p(i) = -imu_omega_B{i}(2,3);
    randwalk_r = randwalk_r + gyro_randwalk_enable * 0.1*(1/t(end))*randn;
    randwalk_q = randwalk_q + gyro_randwalk_enable * 0.1*(1/t(end))*randn;
    randwalk_p = randwalk_p + gyro_randwalk_enable * 0.1*(1/t(end))*randn;
    meas_r(i) = imu_r(i) + sigma_pqr*randn + randwalk_r;
    meas_q(i) = imu_q(i) + sigma_pqr*randn + randwalk_q;
    meas_p(i) = imu_p(i) + sigma_pqr*randn + randwalk_p;
end

figure;
plot(t,meas_r,t,meas_q,t,meas_p,t,imu_r,t,imu_q,t,imu_p);
legend('meas r','meas q','meas p','imu r','imu q','imu p');

% Calculate (body) accelerations and generate noise-corrupted measurements
imu_f_B = {};
for i=1:length(t)
  imu_f_B{i} = imu_R_SB{i}' * [0;0;-g];  %NED, upwards measured acceleration
end

imu_x = NaN*ones(size(t));
imu_y = NaN*ones(size(t));
imu_z = NaN*ones(size(t));
meas_x = NaN*ones(size(t));
meas_y = NaN*ones(size(t));
meas_z = NaN*ones(size(t));
for i=1:length(t)
    imu_x(i) = imu_f_B{i}(1);
    imu_y(i) = imu_f_B{i}(2);
    imu_z(i) = imu_f_B{i}(3);
    meas_x(i) = imu_x(i) + sigma_xyz*randn;
    meas_y(i) = imu_y(i) + sigma_xyz*randn;
    meas_z(i) = imu_z(i) + sigma_xyz*randn;
end

figure;
plot(t,meas_x,t,meas_y,t,meas_z,t,imu_x,t,imu_y,t,imu_z);
legend('meas x','meas y','meas z','imu x','imu y','imu z');

% Calculate heading and generate noise-corrupted measurements
mag_B = {};
for i=1:length(t)
    mag_B{i} = imu_R_SB{i}' * [mag_G;0;0];  %NED, assumed to be true North-aligned and with zero dip at the operating location
end

mag_x = NaN*ones(size(t));
mag_y = NaN*ones(size(t));
mag_z = NaN*ones(size(t));
meas_mag_x = NaN*ones(size(t));
meas_mag_y = NaN*ones(size(t));
meas_mag_z = NaN*ones(size(t));
for i=1:length(t)
    mag_x(i) = mag_B{i}(1);
    mag_y(i) = mag_B{i}(2);
    mag_z(i) = mag_B{i}(3);
    meas_mag_x(i) = mag_x(i) + sigma_mag*randn;
    meas_mag_y(i) = mag_y(i) + sigma_mag*randn;
    meas_mag_z(i) = mag_z(i) + sigma_mag*randn;
end

figure;
plot(t,meas_mag_x,t,meas_mag_y,t,meas_mag_z,t,mag_x,t,mag_y,t,mag_z);
legend('meas mag x','meas mag y','meas mag z','mag x','mag y','mag z');

%% Case 1: Naive "bad" estimates using only accelerometer (& magnetometer) measurements (https://wiki.dfrobot.com/How_to_Use_a_Three-Axis_Accelerometer_for_Tilt_Sensing) %%
bad_accelonly_roll = -atan2(meas_y,-meas_z);
bad_accelonly_pitch = atan2(meas_x,sqrt(meas_y.^2+meas_z.^2));
bad_accelonly_yaw = atan2(-meas_mag_y,meas_mag_x);
bad_accelonly_R_SB = {};
for i=1:length(t)
  bad_accelonly_Eul = [bad_accelonly_yaw(i) bad_accelonly_pitch(i) bad_accelonly_roll(i)];
  bad_accelonly_R_SB{i} = eul2rotm(bad_accelonly_Eul,'ZYX');
end

% For visual comparison
real_yaw = NaN*ones(size(t));
real_pitch = NaN*ones(size(t));
real_roll = NaN*ones(size(t));
for i=1:length(t)
    real_Eul = rotm2eul(imu_R_SB{i},'ZYX');
    real_yaw(i) = real_Eul(1);
    real_pitch(i) = real_Eul(2);
    real_roll(i) = real_Eul(3);
end

figure;
plot(t,bad_accelonly_yaw,t,bad_accelonly_pitch,t,bad_accelonly_roll,t,real_yaw,t,real_pitch,t,real_roll);
legend('bad mag-only yaw','bad accel-only pitch','bad accel-only roll','real yaw','real pitch','real roll');


%% Case 2: Gyro dynamical model prediction (https://www.ethaneade.com/latex2html/lie/node11.html) %%
gyromodel_Sigma = eye(3) * (sigma_pqr*Ts)^2;

gyromodel_R_SB_mu = {};
gyromodel_R_SB_Sigma = {};

gyromodel_R_SB_mu{1} = imu_R_SB{1};  %initial mean
gyromodel_R_SB_Sigma{1} = eye(3) * deg2rad(1e-9);  %initial covariance

for i=2:length(t)
  gyromodel_meas_skewsymmetric = [0          -meas_r(i) meas_q(i);
                                  meas_r(i)  0          -meas_p(i);
                                  -meas_q(i) meas_p(i)  0];
  
  gyromodel_R_SB_mu{i} = gyromodel_R_SB_mu{i-1} * expm(gyromodel_meas_skewsymmetric * Ts);
  gyromodel_R_SB_Sigma{i} = gyromodel_R_SB_Sigma{i-1} + gyromodel_R_SB_mu{i-1} * gyromodel_Sigma * gyromodel_R_SB_mu{i-1}';

  % Alternatively:
  % gyromodel_R_SB_mu_dot = gyromodel_R_SB_mu{i-1} * gyromodel_meas_skewsymmetric;
  % gyromodel_R_SB_mu{i} = gyromodel_R_SB_mu{i-1} + gyromodel_R_SB_mu_dot * Ts;
end

% Calculate "naive" predictions
bad_gyroonly_R_SB = gyromodel_R_SB_mu;
bad_gyroonly_yaw = NaN*ones(size(t));
bad_gyroonly_pitch = NaN*ones(size(t));
bad_gyroonly_roll = NaN*ones(size(t));
for i=1:length(t)
  bad_gyroonly_Eul = rotm2eul(bad_gyroonly_R_SB{i},'ZYX');
  bad_gyroonly_yaw(i) = bad_gyroonly_Eul(1);
  bad_gyroonly_pitch(i) = bad_gyroonly_Eul(2);
  bad_gyroonly_roll(i) = bad_gyroonly_Eul(3);
end

figure;
plot(t,bad_gyroonly_yaw,t,bad_gyroonly_pitch,t,bad_gyroonly_roll,t,real_yaw,t,real_pitch,t,real_roll);
legend('bad gyroonly yaw','bad gyroonly pitch','bad gyroonly roll','real yaw','real pitch','real roll');


%% Case 3: Bayesian Estimation (https://www.ethaneade.com/latex2html/lie/node9.html) %%
gyromodel_Sigma = eye(3) * (sigma_pqr*Ts)^2;
accelobservationmodel_Sigma = eye(3) * deg2rad(10)^2;  %assume direct observation model for rotation with 10 deg sigma on each axis

estim_bayes_R_SB_mu = {};
estim_bayes_R_SB_Sigma = {};

estim_bayes_R_SB_mu{1} = imu_R_SB{1};  %initial mean
estim_bayes_R_SB_Sigma{1} = eye(3) * deg2rad(1e-9);  %initial covariance

for i=2:length(t)
  % Predict
  gyromodel_meas_skewsymmetric = [0          -meas_r(i) meas_q(i);
                                  meas_r(i)  0          -meas_p(i);
                                  -meas_q(i) meas_p(i)  0];
  
  predict_R_SB_mu = estim_bayes_R_SB_mu{i-1} * expm(gyromodel_meas_skewsymmetric * Ts);
  predict_R_SB_Sigma = estim_bayes_R_SB_Sigma{i-1} + estim_bayes_R_SB_mu{i-1} * gyromodel_Sigma * estim_bayes_R_SB_mu{i-1}';

  % Update
  %Sigma_combined = pinv(pinv(predict_R_SB_Sigma) + pinv(accelobservationmodel_Sigma))
  Sigma_combined = predict_R_SB_Sigma - predict_R_SB_Sigma * pinv(predict_R_SB_Sigma + accelobservationmodel_Sigma) * predict_R_SB_Sigma;
  innov = logm(bad_accelonly_R_SB{i} * predict_R_SB_mu');
  innov(1,1) = 0; innov(2,2) = 0; innov(3,3) = 0;  innov(2,1) = -innov(1,2); innov(3,1) = - innov(1,3); innov(3,2) = - innov(2,3); %force skew-symmetric 
  R_SB_combined = expm(Sigma_combined * pinv(accelobservationmodel_Sigma) * innov) * predict_R_SB_mu;

  estim_bayes_R_SB_mu{i} = R_SB_combined;
  estim_bayes_R_SB_Sigma{i} = Sigma_combined;
end

% Calculate Bayesian-fused estimates
bayes_R_SB = estim_bayes_R_SB_mu;
bayes_yaw = NaN*ones(size(t));
bayes_pitch = NaN*ones(size(t));
bayes_roll = NaN*ones(size(t));
for i=1:length(t)
  bayes_Eul = rotm2eul(bayes_R_SB{i},'ZYX');
  bayes_yaw(i) = bayes_Eul(1);
  bayes_pitch(i) = bayes_Eul(2);
  bayes_roll(i) = bayes_Eul(3);
end

figure;
plot(t,bayes_yaw,t,bayes_pitch,t,bayes_roll,t,real_yaw,t,real_pitch,t,real_roll);
legend('bayes yaw','bayes pitch','bayes roll','real yaw','real pitch','real roll');


%% Case 4: Quaternion Extended Kalman Filtering (https://ahrs.readthedocs.io/en/latest/filters/ekf.html & https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles) %%
gyromodel_Sigma = eye(3) * (sigma_pqr)^2;
accelmodel_Sigma = eye(3) * (sigma_xyz)^2;
magmodel_Sigma = eye(3) * (sigma_mag)^2;

estim_kalman_q_SB_mu = {};
estim_kalman_q_SB_Sigma = {};

estim_kalman_q_SB_mu{1} = rotm2quat( imu_R_SB{1} )';  %initial mean
estim_kalman_q_SB_Sigma{1} = eye(4) * deg2rad(1e-9);  %initial covariance

for i=2:length(t)
  % Predict
  qw = estim_kalman_q_SB_mu{i-1}(1);
  qx = estim_kalman_q_SB_mu{i-1}(2);
  qy = estim_kalman_q_SB_mu{i-1}(3);
  qz = estim_kalman_q_SB_mu{i-1}(4);

  gyromodel_meas = [meas_p(i);
                    meas_q(i);
                    meas_r(i)];
  gyromodel_meas_skewsymmetric = [0          -meas_r(i) meas_q(i);
                                  meas_r(i)  0          -meas_p(i);
                                  -meas_q(i) meas_p(i)  0];
  Omega_quat_skewsymmetric = [0               -gyromodel_meas';
                              [gyromodel_meas -gyromodel_meas_skewsymmetric]];

  predict_q_SB_mu = (eye(4) + (Ts/2)*Omega_quat_skewsymmetric) * estim_kalman_q_SB_mu{i-1};

  F_q = (eye(4) + (Ts/2)*Omega_quat_skewsymmetric);
  F_gyromodel_meas = (Ts/2)*[-qx -qy -qz;
                             qw  qz  -qy;
                             -qz qw  qx;
                             qy  -qx qw];
  predict_q_SB_Sigma = F_q * estim_kalman_q_SB_Sigma{i-1} * F_q' + F_gyromodel_meas * gyromodel_Sigma * F_gyromodel_meas';

  % Update
  gx = 0;
  gy = 0;
  gz = -1;  % normalized gravity vector

  rx = 1;  % normalized magnetic North vector
  ry = 0;
  rz = 0;

  qw = predict_q_SB_mu(1);
  qx = predict_q_SB_mu(2);
  qy = predict_q_SB_mu(3);
  qz = predict_q_SB_mu(4);
  predict_R_SB_mu = quat2rotm(predict_q_SB_mu');

  accel_meas = [meas_x(i);
                meas_y(i);
                meas_z(i)];
  accel_meas = accel_meas/norm(accel_meas);  %normalize accel measurement
  predict_accel_meas = predict_R_SB_mu' * [gx;gy;gz];
  innov_accel = accel_meas - predict_accel_meas;
  
  H_accel_q = 2 * [gy*qz-gz*qy   gy*qy+gz*qz          -2*gx*qy+gy*qx-gz*qw  -2*gx*qz+gy*qw+gz*qx;
                   -gx*qz+gz*qx  gx*qy-2*gy*qx+gz*qw  gx*qx+gz*qz           -gx*qw-2*gy*qz+gz*qy;
                   gx*qy-gy*qx   gx*qz-gy*qw-2*gz*qx  gx*qw+gy*qz-2*gz*qy   gx*qx+gy*qy];

  mag_meas = [meas_mag_x(i);
              meas_mag_y(i);
              meas_mag_z(i)];
  mag_meas = mag_meas/norm(mag_meas);  %normalize mag measurement
  predict_mag_meas = predict_R_SB_mu' * [rx;ry;rz];
  innov_mag = mag_meas - predict_mag_meas;

  H_mag_q = 2 * [rx*qw+ry*qz-rz*qy   rx*qx+ry*qy+rz*qz  -rx*qy+ry*qx-rz*qw  -rx*qz+ry*qw+rz*qx;
                 -rx*qz+ry*qw+rz*qx  rx*qy-ry*qx+rz*qw  rx*qx+ry*qy+rz*qz   -rx*qw-ry*qz+rz*qy;
                 rx*qy-ry*qx+rz*qw   rx*qz-ry*qw-rz*qx  rx*qw+ry*qz-rz*qy   rx*qx+ry*qy+rz*qz];

  H_q = [H_accel_q;H_mag_q];
  innov = [innov_accel;innov_mag];
  R = blkdiag(accelmodel_Sigma,magmodel_Sigma);

  S = H_q * predict_q_SB_Sigma * H_q' + R;
  K = predict_q_SB_Sigma * H_q' * pinv(S);
  estim_kalman_q_SB_mu{i} = predict_q_SB_mu + K * innov;
  estim_kalman_q_SB_Sigma{i} = (eye(4) - K * H_q) * predict_q_SB_Sigma;

  estim_kalman_q_SB_mu{i} = estim_kalman_q_SB_mu{i}/norm(estim_kalman_q_SB_mu{i});  %fix to retain normalization
end

% Calculate Kalman-fused estimates
kalman_q_SB = estim_kalman_q_SB_mu;
kalman_yaw = NaN*ones(size(t));
kalman_pitch = NaN*ones(size(t));
kalman_roll = NaN*ones(size(t));
for i=1:length(t)
  kalman_Eul = quat2eul(kalman_q_SB{i}','ZYX');
  kalman_yaw(i) = kalman_Eul(1);
  kalman_pitch(i) = kalman_Eul(2);
  kalman_roll(i) = kalman_Eul(3);
end

% For visual comparison
real_yaw = NaN*ones(size(t));
real_pitch = NaN*ones(size(t));
real_roll = NaN*ones(size(t));
for i=1:length(t)
    real_Eul = rotm2eul(imu_R_SB{i},'ZYX');
    real_yaw(i) = real_Eul(1);
    real_pitch(i) = real_Eul(2);
    real_roll(i) = real_Eul(3);
end

figure;
plot(t,kalman_yaw,t,kalman_pitch,t,kalman_roll,t,real_yaw,t,real_pitch,t,real_roll);
legend('kalman yaw','kalman pitch','kalman roll','real yaw','real pitch','real roll');

%% Plotting
scenario = 1;  % 1=accelonly , 2=gyroonly , 3=bayesian , 4=kalman

figure;
ax = axes;
draw_stepsize = 1;
for i=1:draw_stepsize:length(t)
  cla(ax);
  hold(ax,'on');

  plotTransforms([0 0 0],rotm2quat(imu_R_SB{i}));

%   acc_x_vec = imu_R_SB{i}*[imu_x(i);0;0];
%   quiver3([0],[0],[0], 0.1*[acc_x_vec(1)],0.1*[acc_x_vec(2)],0.1*[acc_x_vec(3)], 'LineWidth',2,'Color',[0.64,0.08,0.18]);  %scale magnitude for visualization
%   acc_y_vec = imu_R_SB{i}*[0;imu_y(i);0];
%   quiver3([0],[0],[0], 0.1*[acc_y_vec(1)],0.1*[acc_y_vec(2)],0.1*[acc_y_vec(3)], 'LineWidth',2,'Color',[0.47,0.67,0.19]);  %scale magnitude for visualization
%   acc_z_vec = imu_R_SB{i}*[0;0;imu_z(i)];
%   quiver3([0],[0],[0], 0.1*[acc_z_vec(1)],0.1*[acc_z_vec(2)],0.1*[acc_z_vec(3)], 'LineWidth',2,'Color',[0.30,0.75,0.93]);  %scale magnitude for visualization
%   acc_vec = acc_x_vec+acc_y_vec+acc_z_vec;
%   quiver3([0],[0],[0], 0.1*[acc_vec(1)],0.1*[acc_vec(2)],0.1*[acc_vec(3)], 'LineWidth',2,'Color',[0.49,0.18,0.56]);  %scale magnitude for visualization

  if scenario == 1
      plotTransforms([0 0 0],rotm2quat(bad_accelonly_R_SB{i}));
  elseif scenario == 2
      plotTransforms([0 0 0],rotm2quat(bad_gyroonly_R_SB{i}));
  elseif scenario == 3
      plotTransforms([0 0 0],rotm2quat(bayes_R_SB{i})); 
  elseif scenario == 4
      plotTransforms([0 0 0],estim_kalman_q_SB_mu{i}');
  end

  hold(ax,'off');
  set(ax,'zdir','reverse','ydir','reverse');  %NED
  set(ax,'xlim',[-1,1],'ylim',[-1,1],'zlim',[-1,1]);
  grid(ax,'on');

  drawnow;

  if i==1; pause; end
end