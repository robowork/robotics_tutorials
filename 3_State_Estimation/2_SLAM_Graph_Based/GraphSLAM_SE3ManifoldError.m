clear all; clc; bdclose all;
load('data/sample_unoptimized_EKF_SLAM_data.mat');

num_poses = length(ksi_estimates);
num_landmarks = length(m_groundtruth)/2; 
num_inputs = length(u_inputs);

colors = {'red','green','blue','magenta','yellow','cyan'};

Dt = 0.2; %static,if changing then for every pose-pose do: ksi_estimates{1 + t}.time - ksi_estimates{1 + t-1}.time;

sigma_v = 0.25; % forward velocity
sigma_s = 0.025; % account for some lateral wheel slippage velocity
sigma_w = deg2rad(2.5); %rotational velocity
Q = blkdiag(sigma_v^2 * Dt, sigma_s^2 * Dt, 1e-12 , 1e-12, 1e-12, sigma_w^2 * Dt);
Q_inv = inv(Q);

sigma_x = 0.5; % landmark x localization %static,if changing then for every pose-landmark estimate based on range
sigma_y = 0.5; % landmark y localization %static,if changing then for every pose-landmark estimate based on range
R = blkdiag(sigma_x^2, sigma_y^2, 1e-12);
R_inv = inv(R);

%%% Align the groundtruth and estimates initial frames of reference %%%
% Groundtruth initial transformation
groundtruth_trans_0 = [ksi_actual{1}.ksi_groundtruth(1);...
                       ksi_actual{1}.ksi_groundtruth(2);
                       0];
groundtruth_rot_0_rotm = axang2rotm([0 0 1 ksi_actual{1}.ksi_groundtruth(3)]);
groundtruth_0_hom = [groundtruth_rot_0_rotm, groundtruth_trans_0;
                     0, 0, 0,                1];
% Estimate initial transformation
estim_trans_0 = [ksi_estimates{1}.ksi_estim(1);...
                 ksi_estimates{1}.ksi_estim(2);
                 0];
estim_rot_0_rotm = axang2rotm([0 0 1 ksi_estimates{1}.ksi_estim(3)]);
estim_0_hom = [estim_rot_0_rotm, estim_trans_0;
               0, 0, 0,          1];
% Bring Estimate onto Groundtruth initial transformation
for tt=0:num_poses-1
    ksi_estim_trans = [ksi_estimates{1 + tt}.ksi_estim(1);...
                       ksi_estimates{1 + tt}.ksi_estim(2);...
                       0];
    ksi_estim_rot_rotm = axang2rotm([0 0 1 ksi_estimates{1 + tt}.ksi_estim(3)]);
    ksi_estim_hom = [ksi_estim_rot_rotm, ksi_estim_trans;
                     0, 0, 0,       1];           
    ksi_estim_hom = groundtruth_0_hom * inv(estim_0_hom) * ksi_estim_hom;
    ksi_estim_trans = ksi_estim_hom(1:3,4);
    ksi_estim_rot_rotm = ksi_estim_hom(1:3,1:3);
    ksi_axang_estim = rotm2axang(ksi_estim_rot_rotm);
    if (ksi_axang_estim(3) < 0); ksi_axang_estim = -ksi_axang_estim; end
    ksi_estimates{1 + tt}.ksi_estim = [ksi_estim_trans(1:2);ksi_axang_estim(4)];
end 

% Unwrap to get continuous yaw numbers (we won't do actual SE(2) math, for this sample to remain simple we'll do "fake-flat-like-space" optimization)
for tt=0:num_poses-1 
    if tt >= 1 
      while ksi_estimates{1 + tt}.ksi_estim(3) - ksi_estimates{1 + tt-1}.ksi_estim(3) > pi
        ksi_estimates{1 + tt}.ksi_estim(3) = ksi_estimates{1 + tt}.ksi_estim(3) - 2*pi;
      end
      while ksi_estimates{1 + tt}.ksi_estim(3) - ksi_estimates{1 + tt-1}.ksi_estim(3) < -pi
        ksi_estimates{1 + tt}.ksi_estim(3) = ksi_estimates{1 + tt}.ksi_estim(3) + 2*pi;
      end
    end
end

h_fig = figure;
axis equal;
set(gca,'xlim',[-10.0 10.0],'ylim',[-10 10],'zlim',[-10 10]);
view(0,90);

hold on; 
for nn=0:num_landmarks - 1
    plot(gca,m_groundtruth(1 + 2*nn),m_groundtruth(1 + 2*nn+1),'Marker','o','Color','black','MarkerSize',10);
    %line(gca,[ksi_groundtruth(1) m_groundtruth(1 + 2*nn)],[ksi_groundtruth(2) m_groundtruth(1 + 2*nn+1)],'LineStyle','--','Color',[0.75,0.75,0.75]);
end
for tt=0:num_poses-1
    time = m_estimates{1 + tt}.time;
    m_estim = m_estimates{1 + tt}.m_estim;
    for nn=0:num_landmarks - 1
        plot(gca,m_estim(1 + 2*nn),m_estim(1 + 2*nn+1),'Marker','*','Color',colors{1 + mod(tt,length(colors))}); 
    end
end
for tt=0:num_poses-1
    time = ksi_actual{1 + tt}.time;
    ksi_groundtruth = ksi_actual{1 + tt}.ksi_groundtruth;
    trans_groundtruth = [ksi_groundtruth(1) ksi_groundtruth(2) 0];
    rot_groundtruth = axang2quat([0 0 1 ksi_groundtruth(3)]);
    h_trans_groundtruth = plotTransforms(trans_groundtruth,rot_groundtruth); 
    plot(gca,ksi_groundtruth(1),ksi_groundtruth(2),'Marker','x','Color',colors{1 + mod(tt,length(colors))}); 
    view(0,90);
end
for tt=0:num_inputs-1
    time = u_inputs{1 + tt}.time;
    u_t = u_inputs{1 + tt}.u_t;
end
for tt=0:num_poses-1
    time = ksi_estimates{1 + tt}.time;
    ksi_estim = ksi_estimates{1 + tt}.ksi_estim;
    z_mi = z_measurements{1 + tt}.z_mi;
    plot(gca,ksi_estim(1),ksi_estim(2),'Marker','o','Color',colors{1 + mod(tt,length(colors))}); 
    line(gca,[ksi_estim(1) ksi_estim(1)+1.0*cos(ksi_estim(3))],[ksi_estim(2) ksi_estim(2)+1.0*sin(ksi_estim(3))],'LineStyle','--','Color','black'); 
    for kk=0:length(z_mi)/2 - 1
      %line(gca,[ksi_estim(1) ksi_estim(1)+z_mi(1 + 2*kk)*cos(ksi_estim(3)+z_mi(1 + 2*kk+1))],[ksi_estim(2) ksi_estim(2)+z_mi(1 + 2*kk)*sin(ksi_estim(3)+z_mi(1 + 2*kk+1))],'LineStyle','--','Color',colors{1 + mod(tt,length(colors))});
    end
end
hold off;

%%% Graph-SLAM %%%
DX = NaN * zeros(6*num_poses+3*num_landmarks , 1); % [t,w](in R^6)*x_t + [t_m](in R^3)*m_n
DX_poses_index = 1;
DX_landmarks_index = 1+6*num_poses;

% List of nodes with SE(3) pose estimates from last Kalman estimate at that time
X = {};
for tt=0:num_poses-1
    X{1 + tt}.T = [ksi_estimates{1 + tt}.ksi_estim(1:2); 0];
    X{1 + tt}.R = axang2rotm([0 0 1 ksi_estimates{1 + tt}.ksi_estim(3)]);
end
% List of nodes with landmarks' locations from final Kalman estimate in the recorded sequence
M = {};
for nn=0:num_landmarks-1
    M{1 + nn}.T = [m_estimates{1 + nn}.m_estim(1:2); 0];
end

% ITERATIVELY
for ITERATION=1:100
        
    H = zeros(6*num_poses+3*num_landmarks , 6*num_poses+3*num_landmarks); 
    b = zeros(6*num_poses+3*num_landmarks , 1); 
    
    % For every pose recorded in our graph
    for t=0:num_poses-1
      
      % Graph-SLAM (28)    
      if t == 0
        % 1-st pose of our graph, add anchoring constraint
        H(DX_poses_index : DX_poses_index +6-1, DX_poses_index : DX_poses_index + 6-1) = 1 * eye(6);
        b(DX_poses_index : DX_poses_index +6-1) = zeros(6,1); 
        
      elseif t >= 1
        u_t = u_inputs{1 + t}.u_t; % t-th control input (corresponds to t-1 applied robot control)
        % i: (t-1)-th pose, j: (t)-th pose
        % error term
        ksi_estim_tm1.T = X{1 + t-1}.T;
        ksi_estim_tm1.R = X{1 + t-1}.R;
        ksi_estim_t.T = X{1 + t}.T;
        ksi_estim_t.R = X{1 + t}.R;
        % (t-1)->(t)-th relative transformation "virtual measurement" from motion model
        z_ij.T = [u_t(1)*Dt; ...
                  0; ...
                  0];
        z_ij.R = axang2rotm([0 0 1 u_t(2)*Dt]);  
        
        % homogeneous transformations
        Ksi_i = [ksi_estim_tm1.R  ksi_estim_tm1.T;...
                 0 0 0                          1];
        Ksi_i_inv = inv(Ksi_i);
        Ksi_j = [ksi_estim_t.R  ksi_estim_t.T;...
                 0 0 0                      1];
        Ksi_i_inv_dot_Ksi_j = Ksi_i_inv * Ksi_j;
        Z_ij = [z_ij.R  z_ij.T;...
                0 0 0        1];
        Z_ij_inv = inv(Z_ij);
        
        % error SE(3) formulation (Z_ij^-1 * (X_i^-1 * X_j))
        E_ij = Z_ij_inv * Ksi_i_inv_dot_Ksi_j;
        e_ij.R = E_ij(1:3,1:3);
        e_ij.T = E_ij(1:3,4);
        
        % https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.10.3.11
        cos_theta = (trace(e_ij.R)-1)/2;
        if cos_theta > 0.999999
            theta = 0;
            d_log_ER_d_ER = [0   0    0 ,    0 0 0.5 ,   0 -0.5 0;...
                             0   0 -0.5 ,    0 0   0 , 0.5    0 0;...
                             0 0.5    0 , -0.5 0   0 ,   0    0 0];
        else
            theta = acos(cos_theta);
            sin_theta = sqrt(1-cos_theta^2);
            alpha = [e_ij.R(3,2)-e_ij.R(2,3);...
                     e_ij.R(1,3)-e_ij.R(3,1);...
                     e_ij.R(2,1)-e_ij.R(1,2)] * (theta*cos_theta-sin_theta)/(4*sin_theta^3);
            beta = theta/(2*sin_theta);
            d_log_ER_d_ER = [alpha(1)    0     0 ,     0 alpha(1) beta ,    0 -beta alpha(1);...
                             alpha(2)    0 -beta ,     0 alpha(2)    0 , beta     0 alpha(2);...
                             alpha(3) beta     0 , -beta alpha(3)    0 ,    0     0 alpha(3)];
        end
        
        % https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.10.3.35
        d_log_E_d_E = [zeros(3,9)    , eye(3);...
                       d_log_ER_d_ER , zeros(3,3)];  
        
        d_Zinv_dot_XiinvXj_d_Zinv = kron(Ksi_i_inv_dot_Ksi_j,eye(3));
        
        % jacobian in eps1
        d_Zinv_dot_expeps1_d_eps1 = [zeros(3,3)        ,  zeros(3,1)      , -Z_ij_inv(1:3,3) ,  Z_ij_inv(1:3,2);...
                                     zeros(3,3)        ,  Z_ij_inv(1:3,3) ,  zeros(3,1)      , -Z_ij_inv(1:3,1);... 
                                     zeros(3,3)        , -Z_ij_inv(1:3,2) ,  Z_ij_inv(1:3,1) ,  zeros(3,1);...
                                     Z_ij_inv(1:3,1:3) ,                     zeros(3,3)                   ];
        
        J_ij_eps1 = d_log_E_d_E * d_Zinv_dot_XiinvXj_d_Zinv * (-d_Zinv_dot_expeps1_d_eps1);
        
        % jacobian in eps2
        d_ZinvXiinvXj_dot_expeps2_d_eps2 = [zeros(3,3)    ,  zeros(3,1)  , -E_ij(1:3,3) ,  E_ij(1:3,2);...
                                            zeros(3,3)    ,  E_ij(1:3,3) ,  zeros(3,1)  , -E_ij(1:3,1);... 
                                            zeros(3,3)    , -E_ij(1:3,2) ,  E_ij(1:3,1) ,  zeros(3,1);...
                                            E_ij(1:3,1:3) ,                 zeros(3,3)               ];

        J_ij_eps2 = d_log_E_d_E * (+d_ZinvXiinvXj_dot_expeps2_d_eps2);

                
        % information matrix and vector
        %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.21
        cos_theta = (trace(e_ij.R)-1)/2;
        if cos_theta > 0.999999
          omega_ER = zeros(3,1);
          %%% V_E_ij_R_inv = eye(3);
        else
          theta = acos(cos_theta);
          sin_theta = sqrt(1-cos_theta^2);
          beta = theta/(2*sin_theta);
          log_ER = beta * (e_ij.R-e_ij.R');  %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.16
          omega_skew = log_ER;  %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.16
          omega_ER = [omega_skew(3,2);omega_skew(1,3);omega_skew(2,1)];  %vee of skew-symmetric
          %%% V_E_ij_R_inv = eye(3) - (1/2)*(omega_skew) + ((1-((theta*cos(theta/2))/(2*sin(theta/2))))/(theta^2))*(omega_skew^2);  %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.26
        end
        %%% t_E_ij = V_E_ij_R_inv * e_ij.T;  %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.24
        %%% log_E_ij = [t_E_ij;omega_E_ij_R];
        log_E_ij = [e_ij.T;omega_ER];  %pseudo-log: %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.29
      
        
        J_ij_big = [J_ij_eps1 , J_ij_eps2];
        H_ij_big = J_ij_big' * Q_inv * J_ij_big;
        b_ij_big = J_ij_big' * Q_inv * log_E_ij;
              
                 
        % update big information matrix and vector
        H(DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1, DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1) = ...
          H(DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1, DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1) + ... 
            H_ij_big(1:6,1:6);
        H(DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) = ...
          H(DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) + ... 
            H_ij_big(1:6,7:12);
        H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1) = ...
          H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1) + ... 
            H_ij_big(7:12,1:6);
        H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) = ...
          H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) + ... 
            H_ij_big(7:12,7:12);
        
        b(DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1) = ...
          b(DX_poses_index + (t-1)*6 : DX_poses_index + (t)*6-1) + ...
            b_ij_big(1:6);
        b(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) = ...
          b(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) + ...
            b_ij_big(7:12);
      end

      % for every landmark of our graph
      z_mi = z_measurements{1 + t}.z_mi; % actual measurements z for t-th robot pose -> landmark_i:1...n
      for n=0:num_landmarks-1
          ra_t = z_mi(1 + n*2 : 1 + (n+1)*2-1); % actual range & bearing measurement
          % i: (t)-th pose, j: (n)-th landmark
          % error term
          ksi_estim_t.T = X{1 + t}.T;
          ksi_estim_t.R = X{1 + t}.R;
          mn_t.T = M{1 + n}.T;
          % (t-1)->(t)-th relative transformation "virtual measurement" from motion model
          z_ij.T = [ra_t(1)*cos(ra_t(2));...
                    ra_t(1)*sin(ra_t(2));...
                    0];
        
          % homogeneous transformations
          Ksi_i = [ksi_estim_t.R  ksi_estim_t.T;...
                   0 0 0                      1];
          Ksi_i_inv = inv(Ksi_i);
          M_j = [mn_t.T;...
                 1];
          Ksi_i_inv_dot_M_j = Ksi_i_inv * M_j;
          Z_ij = [z_ij.T;...
                  1];
                  
          % error SE(3) formulation (Z_ij^-1 * (X_i^-1 * X_j))
          E_ij = Z_ij - Ksi_i_inv_dot_M_j;
           
          % jacobian in eps1
          %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.10.3.21
          dc1 = Ksi_i_inv(1:3,1);
          dc2 = Ksi_i_inv(1:3,2);
          dc3 = Ksi_i_inv(1:3,3);
          dt = Ksi_i_inv(1:3,4);
          dc1_skew = [ 0      -dc1(3)  dc1(2);...
                       dc1(3)  0      -dc1(1);...
                      -dc1(2)  dc1(1)       0];
          dc2_skew = [ 0      -dc2(3)  dc2(2);...
                       dc2(3)  0      -dc2(1);...
                      -dc2(2)  dc2(1)       0];          
          dc3_skew = [ 0      -dc3(3)  dc3(2);...
                       dc3(3)  0      -dc3(1);...
                      -dc3(2)  dc3(1)       0];
          dt_skew = [ 0     -dt(3)  dt(2);...
                      dt(3)  0     -dt(1);...
                     -dt(2)  dt(1)      0];
                  
          d_expeps1Xiiinv_plus_Mj_d_eps1 = kron([mn_t.T' 1],eye(3)) * [zeros(3,3) , -dc1_skew;...
                                                                       zeros(3,3) , -dc2_skew;...
                                                                       zeros(3,3) , -dc3_skew;...
                                                                       eye(3,3)   , -dt_skew];
        
          J_ij_eps1 = -(-d_expeps1Xiiinv_plus_Mj_d_eps1);
        
          % jacobian in eps2
          J_ij_eps2 = -Ksi_i_inv(1:3,1:3);
                
          % information matrix and vector  
          log_E_ij = E_ij(1:3);  %pseudo-log equivalent
          
          J_ij_big = [J_ij_eps1 , J_ij_eps2];
          H_ij_big = J_ij_big' * R_inv * J_ij_big;
          b_ij_big = J_ij_big' * R_inv * log_E_ij;
        
        
          % update big information matrix and vector
          H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) = ...
            H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) + ... 
              H_ij_big(1:6, 1:6);
          H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1) = ...
            H(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1, DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1) + ... 
              H_ij_big(1:6, 7:9);
          H(DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) = ...
            H(DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1, DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) + ... 
              H_ij_big(7:9, 1:6);
          H(DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1, DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1) = ...
            H(DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1, DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1) + ... 
              H_ij_big(7:9, 7:9);
          
          b(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) = ...
            b(DX_poses_index + (t)*6 : DX_poses_index + (t+1)*6-1) + ...
              b_ij_big(1:6);
          b(DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1) = ...
            b(DX_landmarks_index + (n)*3 : DX_landmarks_index + (n+1)*3-1) + ...
              b_ij_big(7:9);
      end  
      
    end

    % Solve the entire sparse linear system to get optimal on-manifold increment
%     Dx_star = - inv(H) * b;

    %%% Option 2: First Marginalize landmarks (Lemma Table 7), then Condition on pose history
    H_xx = H(DX_poses_index : DX_poses_index + 6*num_poses-1 , DX_poses_index : DX_poses_index + 6*num_poses-1);
    H_mm = H(DX_landmarks_index : DX_landmarks_index + 3*num_landmarks-1 , DX_landmarks_index : DX_landmarks_index + 3*num_landmarks-1);
    H_xm = H(DX_poses_index : DX_poses_index + 6*num_poses-1 , DX_landmarks_index : DX_landmarks_index + 3*num_landmarks-1);
    H_mx = H(DX_landmarks_index : DX_landmarks_index + 3*num_landmarks-1 , DX_poses_index : DX_poses_index + 6*num_poses-1);

    b_x = b(DX_poses_index : DX_poses_index + 6*num_poses-1);
    b_m = b(DX_landmarks_index : DX_landmarks_index + 3*num_landmarks-1);

    % Note for more accurate implementation from GraphSLAM (39):
    % The matrix H_mm is block-diagonal. This follows from the way it is constructed, in particular the absence of any links between pairs of features.
    % This makes the inversion efficient: H_mm_inv = Sum_j( J_j' inv(H_jj) J_j)
    H_mm_inv = inv(H_mm);

    % Note for more accurate implementation from GraphSLAM (Table 4): "for EACH feature do ... let τ(j) be the set of all poses xt at which j was observed"
    Marg_H_xx = H_xx - H_xm * H_mm_inv * H_mx;
    Marg_b_x = b_x - H_xm * H_mm_inv * b_m;

    x_poses_star = - inv(Marg_H_xx) * Marg_b_x;

    % Condition landmarks on pose history
    Cond_H_mm = H_mm;
    
    Cond_b_m = b_m + H_mx * x_poses_star;
    
    x_landmarks_star = - inv(Cond_H_mm) * Cond_b_m;
    
    % Stack up poses and landmarks 
    Dx_star = [x_poses_star;x_landmarks_star];
    
    % Update with on-manifold increment to get new best estimates
    for tt=0:num_poses-1 
      dx_star = Dx_star(DX_poses_index + (tt)*6 : DX_poses_index + (tt+1)*6-1);
        
      theta = norm(dx_star(4:6));
      cos_theta = cos(theta);
      sin_theta = sqrt(1-cos_theta^2);
      if cos_theta > 0.999999
        %%% V = eye(3);
        %%% e_VMatrix = [eye(3)  V*dx_star(1:3);
        %%%              0 0 0                1];
        e_VMatrix = [eye(3)  dx_star(1:3);
                     0 0 0              1];  %pseudo-exp: %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.27
      else
        omegaSkew = [0          -dx_star(6)  dx_star(5);
                     dx_star(6)  0          -dx_star(4);
                    -dx_star(5)  dx_star(4)  0];
        e_omegaSkew = eye(3) + ( sin_theta/theta )*(omegaSkew) + ( (1-cos_theta)/(theta^2) )*(omegaSkew^2);
        %%% V = eye(3) + ( (1-cos_theta)/(theta^2) )*(omegaSkew) + ( (theta-sin_theta)/(theta^3) )*(omegaSkew^2);
        %%% e_VMatrix = [e_omegaSkew  V*dx_star(1:3);
        %%%              0 0 0                     1];
        e_VMatrix = [e_omegaSkew  dx_star(1:3);
                     0 0 0                   1];  %pseudo-exp: %https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf#equation.9.4.27
      end
      
      X_hom = [X{1 + tt}.R  X{1 + tt}.T;...
               0 0 0                  1];
      X_hom = X_hom * e_VMatrix;  % Box-plus operation
      X{1 + tt}.R = X_hom(1:3,1:3);
      X{1 + tt}.T = X_hom(1:3,4);
    end
    for nn=0:num_landmarks-1
      dx_star = Dx_star(DX_landmarks_index + (nn)*3 : DX_landmarks_index + (nn+1)*3-1); 
      
      M{1 + nn}.T = M{1 + nn}.T + dx_star;  % this is just flat R^3
    end
    
    % Plot updated optimized estimate every mod number of iterations
    if mod(ITERATION,1) == 0
        pause;
        hold on;
        cla(gca);
            trans_0 = [X{1 + 0}.T(1) - ksi_actual{1}.ksi_groundtruth(1);...
                       X{1 + 0}.T(2) - ksi_actual{1}.ksi_groundtruth(2);...
                       0];
            rot_0_axang_tmp = rotm2axang(X{1 + 0}.R);
            if rot_0_axang_tmp(3)<0; rot_0_axang_tmp = -rot_0_axang_tmp; end
            rot_0_rotm = axang2rotm([0 0 1 rot_0_axang_tmp(4)-ksi_actual{1}.ksi_groundtruth(3)]);
            hom_0 = [rot_0_rotm, trans_0;
                     0, 0, 0,    1];                 
        for nn=0:num_landmarks-1
          plot(gca,m_groundtruth(1 + 2*nn),m_groundtruth(1 + 2*nn+1),'Marker','o','Color','black','MarkerSize',10);
          m_estim = M{1 + nn}.T;
              hom_estim = [eye(3)  m_estim;
                           0 0 0         1];
              hom_estim = inv(hom_0) * hom_estim;
              m_estim = [hom_estim(1,4);hom_estim(2,4)];
          plot(gca,m_estim(1),m_estim(2),'Marker','*','Color','magenta'); 
        end
        for tt=0:num_poses-1
          ksi_groundtruth = ksi_actual{1 + tt}.ksi_groundtruth;
          trans_groundtruth = [ksi_groundtruth(1) ksi_groundtruth(2) 0];
          rot_groundtruth = axang2quat([0 0 1 ksi_groundtruth(3)]);
          h_trans_groundtruth = plotTransforms(trans_groundtruth,rot_groundtruth); 
          plot(gca,ksi_groundtruth(1),ksi_groundtruth(2),'Marker','x','Color',colors{1 + mod(tt,length(colors))}); 
          view(0,90);
              hom_estim = [X{1 + tt}.R, X{1 + tt}.T;
                           0 0 0                  1];
              hom_estim = inv(hom_0) * hom_estim;
              axang_estim = rotm2axang(hom_estim(1:3,1:3));
              if (axang_estim(3) < 0); axang_estim = -axang_estim; end
              ksi_estim = [hom_estim(1,4);hom_estim(2,4);axang_estim(4)];
          plot(gca,ksi_estim(1),ksi_estim(2),'Marker','o','Color',colors{1 + mod(tt,length(colors))}); 
          line(gca,[ksi_estim(1) ksi_estim(1)+1.0*cos(ksi_estim(3))],[ksi_estim(2) ksi_estim(2)+1.0*sin(ksi_estim(3))],'LineStyle','--','Color','black'); 
        end
        hold off;
    end
    
end


