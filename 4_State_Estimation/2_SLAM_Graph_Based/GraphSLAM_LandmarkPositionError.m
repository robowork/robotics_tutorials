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
Q = [sigma_v^2 * Dt,              0,              0;
     0             , sigma_s^2 * Dt,              0;
     0             ,              0, sigma_w^2 * Dt];
Q_inv = inv(Q);

sigma_x = 0.5; % landmark x localization %static,if changing then for every pose-landmark estimate based on range
sigma_y = 0.5; % landmark y localization %static,if changing then for every pose-landmark estimate based on range
R = [sigma_x^2,    0;
     0        ,    sigma_y^2];
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
set(h_fig, 'Color','w');
grid on;
axis equal;
set(gca,'xlim',1.75*[-10.0 10.0],'ylim',1.75*[-10 10],'zlim',1.75*[-10 10]);
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
X = NaN * zeros(3*num_poses+2*num_landmarks , 1); % [x,y,theta]*x_t + [mx,my]*m_n
X_poses_index = 1;
X_landmarks_index = 1+3*num_poses;

% Initialize all nodes with pose estimates from last Kalman estimate at that time
for tt=0:num_poses-1
    ksi_estim = ksi_estimates{1 + tt}.ksi_estim;
    X(X_poses_index + tt*3 : X_poses_index + (tt+1)*3-1) = ksi_estim;
end
% Initialize all nodes with landmarks' locations from final Kalman estimate in the recorded sequence
if num_landmarks >= 1
    tt = num_poses-1;
    m_estim = m_estimates{1 + tt}.m_estim;
    X(X_landmarks_index : X_landmarks_index+2*num_landmarks-1) = m_estim;
end

% ITERATIVELY
for ITERATION=1:1:100
    
    H = zeros(3*num_poses+2*num_landmarks , 3*num_poses+2*num_landmarks); 
    b = zeros(3*num_poses+2*num_landmarks , 1); 
    
    % For every pose recorded in our graph
    for t=0:num_poses-1
      %ksi_estim = X(X_poses_index + t*3 : X_poses_index + (t+1)*3-1);
      %m_estim = X(X_landmarks_index : X_landmarks_index+2*num_landmarks-1);
      
      % Graph-SLAM (28)    
      if t == 0
        % 1-st pose of our graph, add anchoring constraint
        H(X_poses_index : X_poses_index +3-1, X_poses_index : X_poses_index + 3-1) = 1 * eye(3);
        b(X_poses_index : X_poses_index +3-1) = zeros(3,1);  %1 * ksi_estimates{1 + 0}.ksi_estim; 
        
      elseif t >= 1
        u_t = u_inputs{1 + t}.u_t; % t-th control input (corresponds to t-1 applied robot control)
        % i: (t-1)-th pose, j: (t)-th pose
        % error term
        ksi_estim_tm1 = X(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1);
        ksi_estim_t   = X(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1);
        z_ij = [u_t(1)*cos(ksi_estim_tm1(3))*Dt; ...
                u_t(1)*sin(ksi_estim_tm1(3))*Dt; ...
                u_t(2)*Dt]; % (t-1)->(t)-th relative transformation "virtual measurement" from motion model
        e_xy    = -(z_ij(1:2)) + (-ksi_estim_tm1(1:2) + ksi_estim_t(1:2));
        e_theta = -(z_ij(3)) + (-ksi_estim_tm1(3) + ksi_estim_t(3));
        e_ij = [e_xy;...
                e_theta];  % the Euclidian equivalent of SE(3) formulation (Z_ij^-1 * (X_i^-1 * X_j))
            
        % jacobian of error: 
        J_ij = [-1,  0,  u_t(1)*sin(ksi_estim_tm1(3))*Dt, 1, 0, 0;
                 0, -1, -u_t(1)*cos(ksi_estim_tm1(3))*Dt, 0, 1, 0;
                 0,  0, -1                              , 0, 0, 1];

        H_ij = J_ij' * Q_inv * J_ij;
        b_ij = J_ij' * Q_inv * e_ij;
        
        % update big information matrix and vector
        H(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1, X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1) = ...
          H(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1, X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1) + ... 
            H_ij(1:3, 1:3);
        H(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) = ...
          H(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) + ... 
            H_ij(1:3, 4:6);
        H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1) = ...
          H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1) + ... 
            H_ij(4:6, 1:3);
        H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) = ...
          H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) + ... 
            H_ij(4:6, 4:6);
        
        b(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1) = ...
          b(X_poses_index + (t-1)*3 : X_poses_index + (t)*3-1) + ...
            b_ij(1:3);
        b(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) = ...
          b(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) + ...
            b_ij(4:6);
      end

      % for every landmark of our graph
      z_mi = z_measurements{1 + t}.z_mi; % actual measurements z for t-th robot pose -> landmark_i:1...n
      for n=0:num_landmarks-1
          % i: (t)-th pose, j: (n)-th landmark
          % error term
          ksi_estim_t   = X(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1);
          mn_t = X(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1);
          ra_t = z_mi(1 + n*2 : 1 + (n+1)*2-1); % actual range & bearing measurement
          z_ij = [ra_t(1)*cos(ksi_estim_t(3)+ra_t(2)); ...
                  ra_t(1)*sin(ksi_estim_t(3)+ra_t(2));]; % (t)->(n)-th relative transformation "measurement" from sensor model
          e_ij = -(z_ij) + (-ksi_estim_t(1:2) + mn_t);  % the Euclidian equivalent of SE(3) formulation (Z_ij^-1 * (X_i^-1 * X_j))

          % jacobian of error:
          J_ij = [-1,  0,  ra_t(1)*sin(ksi_estim_t(3)+ra_t(2)), 1, 0;
                   0, -1, -ra_t(1)*cos(ksi_estim_t(3)+ra_t(2)), 0, 1];

          H_ij = J_ij' * R_inv * J_ij;
          b_ij = J_ij' * R_inv * e_ij;
    
          % update big information matrix and vector
          H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) = ...
            H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) + ... 
              H_ij(1:3, 1:3);
          H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1) = ...
            H(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1, X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1) + ... 
              H_ij(1:3, 4:5);
          H(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) = ...
            H(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1, X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) + ... 
              H_ij(4:5, 1:3);
          H(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1, X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1) = ...
            H(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1, X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1) + ... 
              H_ij(4:5, 4:5);
          
          b(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) = ...
            b(X_poses_index + (t)*3 : X_poses_index + (t+1)*3-1) + ...
              b_ij(1:3);
          b(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1) = ...
            b(X_landmarks_index + (n)*2 : X_landmarks_index + (n+1)*2-1) + ...
              b_ij(4:5);
      end  
    end

    %%% Option 1: Solve the entire sparse linear system
    % x_star = - inv(H) * b;

    %%% Option 2: First Marginalize landmarks (Lemma Table 7), then Condition on pose history
    H_xx = H(X_poses_index : X_poses_index + 3*num_poses-1 , X_poses_index : X_poses_index + 3*num_poses-1);
    H_mm = H(X_landmarks_index : X_landmarks_index + 2*num_landmarks-1 , X_landmarks_index : X_landmarks_index + 2*num_landmarks-1);
    H_xm = H(X_poses_index : X_poses_index + 3*num_poses-1 , X_landmarks_index : X_landmarks_index + 2*num_landmarks-1);
    H_mx = H(X_landmarks_index : X_landmarks_index + 2*num_landmarks-1 , X_poses_index : X_poses_index + 3*num_poses-1);

    b_x = b(X_poses_index : X_poses_index + 3*num_poses-1);
    b_m = b(X_landmarks_index : X_landmarks_index + 2*num_landmarks-1);

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
    x_star = [x_poses_star;x_landmarks_star];
    
    
    % Update to get new best estimates
    X = X + x_star;   
    % Unwrap to get continuous yaw numbers (we won't do actual SE(2) math, for this sample to remain simple we'll do "fake-flat-like-space" optimization)
    for tt=0:num_poses-1 
        if tt >= 1 
          while X(X_poses_index + 3*tt + 2) - X(X_poses_index + 3*(tt-1) + 2) > pi
            X(X_poses_index + 3*tt + 2) = X(X_poses_index + 3*tt + 2) - 2*pi;
          end
          while X(X_poses_index + 3*tt + 2) - X(X_poses_index + 3*(tt-1) + 2) < -pi
            X(X_poses_index + 3*tt + 2) = X(X_poses_index + 3*tt + 2) + 2*pi;
          end
        end
    end
    

    % Plot updated optimized estimate every mod number of iterations
    if mod(ITERATION,1) == 0
        pause;
        hold on;
        cla(gca);
            trans_0 = [X(X_poses_index + 0) - ksi_actual{1}.ksi_groundtruth(1);...
                       X(X_poses_index + 1) - ksi_actual{1}.ksi_groundtruth(2);...
                       0];
            rot_0 = X(X_poses_index + 2) - ksi_actual{1}.ksi_groundtruth(3);
            rot_0_rotm = axang2rotm([0 0 1 rot_0]);
            hom_0 = [rot_0_rotm, trans_0;
                     0, 0, 0,    1];
        for ii=0:num_landmarks-1
          plot(gca,m_groundtruth(1 + 2*ii),m_groundtruth(1 + 2*ii+1),'Marker','o','Color','black','MarkerSize',10);
          m_estim = X(X_landmarks_index + 2*ii : X_landmarks_index + 2*(ii+1)-1);
              hom_estim = [eye(3),  [m_estim;0];
                           0, 0, 0, 1];
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
          ksi_estim = X(X_poses_index + 3*tt : X_poses_index + 3*(tt+1)-1);
              hom_estim = [axang2rotm([0 0 1 ksi_estim(3)]), [ksi_estim(1:2);0];
                           0, 0, 0,                          1];
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


