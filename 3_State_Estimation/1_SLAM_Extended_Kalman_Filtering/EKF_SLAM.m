global u_t; %Input forward and turn velocity
global z_mi; %Measurments vector for all landmarks
global h_timer; %Handle for timer
global h_fig; %Handle for figure
global h_trans_groundtruth;
global m_groundtruth;
global m_estim;
global diffDriveModel;
global ksi_estim;
global ksi_groundtruth;
global S_ksi_estim;
global S_estim;
global Q;
global R;
global sigma_Rwheel;
global sigma_Lwheel;
global sigma_range;
global sigma_angle;
global dt;
global time;
bdclose all; close all;

m_0_groundtruth = [3;2];  %landmark 0
m_1_groundtruth = [5;-1]; %landmark 1
m_2_groundtruth = [-5;-1]; %landmark 2
m_3_groundtruth = [-4;-5]; %landmark 3
m_4_groundtruth = [-3;5]; %landmark 4
m_groundtruth = [m_0_groundtruth;
                 m_1_groundtruth;
                 m_2_groundtruth;
                 m_3_groundtruth;
                 m_4_groundtruth];  %landmarks

sigma_Rwheel = 1*deg2rad(5);
sigma_Lwheel = 1*deg2rad(5);

sigma_range = 2*0.25;
sigma_angle = 2*deg2rad(2.5);
   
Q = 1*[0.1,0;
     0,deg2rad(1)];
R = 1*[0.25,0;
     0,deg2rad(2.5)];

save('SLAM_data','m_groundtruth', 'sigma_Rwheel', 'sigma_Lwheel', 'sigma_range', 'sigma_angle', 'Q', 'R');

% m_0_estim = 5 * m_0_groundtruth./norm(m_0_groundtruth); %initialize at 5m range
% m_1_estim = 5 * m_1_groundtruth./norm(m_1_groundtruth); %initialize at 5m range
% m_2_estim = 5 * m_2_groundtruth./norm(m_2_groundtruth); %initialize at 5m range
% m_3_estim = 5 * m_3_groundtruth./norm(m_3_groundtruth); %initialize at 5m range
% m_4_estim = 5 * m_4_groundtruth./norm(m_4_groundtruth); %initialize at 5m range
m_0_estim = m_0_groundtruth;
m_1_estim = m_1_groundtruth;
m_2_estim = m_2_groundtruth;
m_3_estim = m_3_groundtruth;
m_4_estim = m_4_groundtruth;
m_estim = [m_0_estim;
           m_1_estim;
           m_2_estim;
           m_3_estim;
           m_4_estim];  %landmarks
       
ksi_0 = [0;0;0]; %initial state [x;y;theta]
ksi_groundtruth = ksi_0;
ksi_estim = ksi_0;
S_ksi_estim = [0.001,0,0;
               0,0.001,0;
               0,0,0.001];
S_m_new = [0.5,0;
           0,0.1];
S_estim = [ S_ksi_estim , zeros(3,length(m_groundtruth));
            zeros(3,length(m_groundtruth))', blkdiag(S_m_new,S_m_new,S_m_new,S_m_new,S_m_new)];

u_t = [0;0];  %Input forward and turn velocity
z_mi = NaN*zeros(length(m_groundtruth) , 1);  %Measurments vector for all landmarks
dt = 0.20; 
time = 0;

diffDriveModel = differentialDriveKinematics("VehicleInputs","VehicleSpeedHeadingRate");  %Vehicle speed and heading angular velocity,specified in meters per second and radians per second respectively. 
diffDriveModel.TrackWidth = 0.5;  %The vehicle track width refers to the distance between the wheels, or the axle length, specified in meters
diffDriveModel.WheelSpeedRange = [-10 10]*2*pi;  %The vehicle speed range is a two-element vector that provides the minimum and maximum vehicle speeds, [MinSpeed MaxSpeed], specified in meters per second.
diffDriveModel.WheelRadius = 0.25;  %The wheel radius of the vehicle, specified in meters.
%[t,pose] = ode45(@(t,pose)derivative(diffDriveModel,y,inputs),tspan,initialState);

h_fig = figure('CloseRequestFcn',@closeFigureReq); 

for i=0:length(m_groundtruth)/2 - 1
    plot(m_groundtruth(1 + 2*i),m_groundtruth(1 + 2*i+1),'o');
end

trans = [ksi_groundtruth(1) ksi_groundtruth(2) 0];
rot = axang2quat([0 0 1 ksi_groundtruth(3)]);
h_trans_groundtruth = plotTransforms(trans,rot); %,'MeshFilePath','groundvehicle.stl'); %,"MeshColor","r");
axis equal;
set(gca,'xlim',[-10.0 10.0],'ylim',[-10 10],'zlim',[-10 10]);
view(0,90);

hold on;
for i=0:length(m_groundtruth)/2 - 1
    plot(m_groundtruth(1 + 2*i),m_groundtruth(1 + 2*i+1),'o');
end
hold off;

set(h_fig,'KeyPressFcn',@keyInput);

h_timer = timer('ExecutionMode','fixedRate','Period',dt,'TimerFcn',@updateDisplay);
start(h_timer);


function updateDisplay(hObject, eventdata)
   global u_t;
   global z_mi;
   global h_fig;
   global h_trans_groundtruth;
   global m_groundtruth;
   global diffDriveModel;
   global ksi_groundtruth;
   global ksi_estim;
   global m_estim;
   global S_ksi_estim;
   global S_estim;
   global Q;
   global R;
   global sigma_Rwheel;
   global sigma_Lwheel;
   global sigma_range;
   global sigma_angle;
   global dt;
   global time;
   
   time = time + dt
   %disp('refresh');
   
   if abs(u_t(1)) > 0; if abs(u_t(1)) <= 0.05; u_t(1) = 0; else; u_t(1) = u_t(1) -sign(u_t(1))*  0.05; end; end
   if abs(u_t(2)) > 0; if abs(u_t(2)) <= 0.01; u_t(2) = 0; else; u_t(2) = u_t(2) -sign(u_t(2))* 0.01; end; end
   
   u_R = u_t(1) + u_t(2) * (diffDriveModel.TrackWidth/2)  +  max([-2*sigma_Rwheel,min([2*sigma_Rwheel,random('Normal',0,sigma_Rwheel)])]);
   u_L = u_t(1) - u_t(2) * (diffDriveModel.TrackWidth/2)  +  max([-2*sigma_Lwheel,min([2*sigma_Lwheel,random('Normal',0,sigma_Lwheel)])]);
   
   u_vel = (u_R + u_L)/2;
   u_rot = (u_R - u_L)/(2 * (diffDriveModel.TrackWidth/2));
   x_dot = u_vel * cos( ksi_groundtruth(3) );
   y_dot = u_vel * sin( ksi_groundtruth(3) );
   theta_dot = u_rot;
   
   ksi_groundtruth(1) = ksi_groundtruth(1) + x_dot * dt;
   ksi_groundtruth(2) = ksi_groundtruth(2) + y_dot * dt;
   ksi_groundtruth(3) = ksi_groundtruth(3) + theta_dot * dt;
   
   trans = [ksi_groundtruth(1) ksi_groundtruth(2) 0];
   rot = axang2quat([0 0 1 ksi_groundtruth(3)]);
   h_trans_groundtruth = plotTransforms(trans,rot);
   axis equal;
   set(gca,'xlim',[-10.0 10.0],'ylim',[-10 10],'zlim',[-10 10]);
   view(0,90);

   hold on;
   for i=0:length(m_groundtruth)/2 - 1
       plot(gca,m_groundtruth(1 + 2*i),m_groundtruth(1 + 2*i+1),'o');
       line(gca,[ksi_groundtruth(1) m_groundtruth(1 + 2*i)],[ksi_groundtruth(2) m_groundtruth(1 + 2*i+1)],'LineStyle','--','Color',[0.75,0.75,0.75]);
   end
   hold off;
   
   z_mi_dx_groundtruth = NaN*zeros(length(m_groundtruth)/2 , 1);
   z_mi_dy_groundtruth = NaN*zeros(length(m_groundtruth)/2 , 1);
   z_mi_dr_groundtruth = NaN*zeros(length(m_groundtruth)/2 , 1);
   z_mi = NaN*zeros(length(m_groundtruth) , 1);
   for i=0:length(m_groundtruth)/2 - 1
       z_mi_dx_groundtruth(1 + i) = m_groundtruth(1 + 2*i)-ksi_groundtruth(1);
       z_mi_dy_groundtruth(1 + i) = m_groundtruth(1 + 2*i+1)-ksi_groundtruth(2);
       z_mi_dr_groundtruth(1 + i) = norm([z_mi_dx_groundtruth(1 + i),z_mi_dy_groundtruth(1 + i)]);
       z_mi(1 + 2*i)   = z_mi_dr_groundtruth(1+i) + max([-2*sigma_range,min([2*sigma_range,random('Normal',0,sigma_range)])]);
       z_mi(1 + 2*i+1) = atan2(z_mi_dy_groundtruth(1 + i),z_mi_dx_groundtruth(1 + i))-ksi_groundtruth(3) + max([-2*sigma_angle,min([2*sigma_angle,random('Normal',0,sigma_angle)])]);
   end
   
   
   ksi_hat = ksi_estim;
   ksi_hat(1) = ksi_hat(1) + dt * u_t(1)*cos(ksi_estim(3) + 0.5*u_t(2)*dt);  %(-u_t(1)/u_t(2))*sin(ksi_estim(3)) + (u_t(1)/u_t(2))*sin(ksi_estim(3) + u_t(2)*dt);
   ksi_hat(2) = ksi_hat(2) + dt * u_t(1)*sin(ksi_estim(3) + 0.5*u_t(2)*dt);  %(u_t(1)/u_t(2))*cos(ksi_estim(3)) + (-u_t(1)/u_t(2))*cos(ksi_estim(3) + u_t(2)*dt);
   ksi_hat(3) = ksi_hat(3) + dt * u_t(2);
   
   F_ksi = [1, 0, dt * -u_t(1)*sin(ksi_estim(3) + 0.5*u_t(2)*dt);
            0, 1, dt * u_t(1)*cos(ksi_estim(3) + 0.5*u_t(2)*dt);
            0, 0, 1];
   F_u = [dt * cos(ksi_estim(3) + 0.5*u_t(2)*dt), dt * u_t(1)*-0.5*dt*sin(ksi_estim(3) + 0.5*u_t(2)*dt);
          dt * sin(ksi_estim(3) + 0.5*u_t(2)*dt), dt * u_t(1)*0.5*dt*cos(ksi_estim(3) + 0.5*u_t(2)*dt);
          0                                     , dt];
   S_ksi_hat = F_ksi * S_estim(1:3,1:3) * F_ksi' + F_u * Q * F_u';
   
   
   z_mi_dx = NaN*zeros(length(m_groundtruth)/2 , 1);
   z_mi_dy = NaN*zeros(length(m_groundtruth)/2 , 1);
   z_mi_dr = NaN*zeros(length(m_groundtruth)/2 , 1);
   z_mi_hat = NaN*zeros(length(m_groundtruth) , 1);
   for i=0:length(m_groundtruth)/2 - 1
       z_mi_dx(1 + i) = m_estim(1 + 2*i)-ksi_hat(1);
       z_mi_dy(1 + i) = m_estim(1 + 2*i+1)-ksi_hat(2);
       z_mi_dr(1 + i) = norm([z_mi_dx(1 + i), z_mi_dy(1 + i)]);
       z_mi_hat(1 + 2*i)   = z_mi_dr(1+i);
       z_mi_hat(1 + 2*i+1) = atan2(z_mi_dy(1+i),z_mi_dx(1+i))-ksi_hat(3);
   end   
   
   H_mi = NaN*zeros(length(m_groundtruth), 5);
   for i=0:length(m_groundtruth)/2 - 1
       H_mi(1 + 2*i : 1 + 2*i+1 , :) = [-z_mi_dx(1+i)/z_mi_dr(1+i)    , -z_mi_dy(1+i)/z_mi_dr(1+i)    , 0 ,  z_mi_dx(1+i)/z_mi_dr(1+i)    , z_mi_dy(1+i)/z_mi_dr(1+i);
                                         z_mi_dy(1+i)/(z_mi_dr(1+i)^2), -z_mi_dx(1+i)/(z_mi_dr(1+i)^2), -1, -z_mi_dy(1+i)/(z_mi_dr(1+i)^2), z_mi_dx(1+i)/(z_mi_dr(1+i)^2)];
   end   
   
   
   % Landmark i
   for i=0:length(m_groundtruth)/2 - 1
       S_hat = zeros(5,5);
       S_hat(1   : 3   , 1   : 3)   = S_ksi_hat;
       S_hat(3+1 : 3+2 , 3+1 : 3+2) = S_estim(3 + (1 + 2*i) : 3 + (1 + 2*i+1) , 3 + (1 + 2*i) : 3 + (1 + 2*i+1));
       S_hat(1   : 3   , 3+1 : 3+2) = S_estim(1             : 3               , 3 + (1 + 2*i) : 3 + (1 + 2*i+1));
       S_hat(3+1 : 3+2 , 1   : 3)   = S_estim(3 + (1 + 2*i) : 3 + (1 + 2*i+1) , 1             : 3);
       K_mi = S_hat * H_mi(1 + 2*i : 1 + 2*i+1 , :)' / (H_mi(1 + 2*i : 1 + 2*i+1 , :) * S_hat * H_mi(1 + 2*i : 1 + 2*i+1 , :)' + R);

       x_hat = [ksi_hat;
                m_estim(1 + 2*i : 1 + 2*i+1)];
       x_upd = x_hat + K_mi * ( z_mi(1 + 2*i : 1 + 2*i+1) - z_mi_hat(1 + 2*i : 1 + 2*i+1) );

       S_upd = ( eye(3+2) - K_mi * H_mi(1 + 2*i : 1 + 2*i+1 , :) ) * S_hat;

       ksi_estim                    = x_upd(1   : 3   , 1);
       m_estim(1 + 2*i : 1 + 2*i+1) = x_upd(3+1 : 3+2 , 1);
       S_estim(1             : 3               , 1             : 3)               = S_upd(1   : 3   , 1   : 3);
       S_estim(3 + (1 + 2*i) : 3 + (1 + 2*i+1) , 3 + (1 + 2*i) : 3 + (1 + 2*i+1)) = S_upd(3+1 : 3+2 , 3+1 : 3+2);
       S_estim(1             : 3               , 3 + (1 + 2*i) : 3 + (1 + 2*i+1)) = S_upd(1   : 3   , 3+1 : 3+2);
       S_estim(3 + (1 + 2*i) : 3 + (1 + 2*i+1) , 1             : 3)               = S_upd(3+1 : 3+2 , 1   : 3);
   end  
   
   % Update globals too
   S_ksi_estim = S_estim(1:3,1:3);
   
   
   % ksi_estim = ksi_hat;
   hold on;
       [eigenvec,eigenval] = eig(S_ksi_estim(1:2,1:2));
       [eigenval_sorted,eigenval_sorted_ind] = sort(diag(eigenval),'descend');
       eigenvec_sorted = eigenvec(:,eigenval_sorted_ind);
       g = hgtransform;
       if eigenval_sorted(1) < 1e-3; eigenval_sorted(1) = 1e-3; end
       if eigenval_sorted(2) < 1e-3; eigenval_sorted(2) = 1e-3; end
       rectangle(gca,'Parent',g,'Position',[-sqrt(5.991*eigenval_sorted(1)) -sqrt(5.991*eigenval_sorted(2)) 2*sqrt(5.991*eigenval_sorted(1)) 2*sqrt(5.991*eigenval_sorted(2))],'Curvature',[1,1],'EdgeColor','magenta');
       g.Matrix=makehgtform('translate',[ksi_estim(1) ksi_estim(2) 0],'zrotate',atan2(eigenvec_sorted(2,1),eigenvec_sorted(1,1)));

       line(gca,[ksi_estim(1) ksi_estim(1)+1.0*cos(ksi_estim(3))],[ksi_estim(2) ksi_estim(2)+1.0*sin(ksi_estim(3))],'Color','red');
         
       mi_g = {};
       for i=0:length(m_groundtruth)/2 - 1
          line(gca,[ksi_estim(1) ksi_estim(1)+z_mi(1 + 2*i)*cos(ksi_estim(3)+z_mi(1 + 2*i+1))],[ksi_estim(2) ksi_estim(2)+z_mi(1 + 2*i)*sin(ksi_estim(3)+z_mi(1 + 2*i+1))],'LineStyle','--','Color','magenta');
                 
          [mi_eigenvec,mi_eigenval] = eig(S_estim(3 + (1 + 2*i) : 3 + (1 + 2*i+1) , 3 + (1 + 2*i) : 3 + (1 + 2*i+1)));
          [mi_eigenval_sorted,mi_eigenval_sorted_ind] = sort(diag(mi_eigenval),'descend');
          mi_eigenvec_sorted = mi_eigenvec(:,mi_eigenval_sorted_ind);
          mi_g{1 + i} = hgtransform;
          if mi_eigenval_sorted(1) < 1e-3; mi_eigenval_sorted(1) = 1e-3; end
          if mi_eigenval_sorted(2) < 1e-3; mi_eigenval_sorted(2) = 1e-3; end
          rectangle(gca,'Parent',mi_g{1 + i},'Position',[-sqrt(5.991*mi_eigenval_sorted(1)) -sqrt(5.991*mi_eigenval_sorted(2)) 2*sqrt(5.991*mi_eigenval_sorted(1)) 2*sqrt(5.991*mi_eigenval_sorted(2))],'Curvature',[1,1],'EdgeColor','black');
          mi_g{1 + i}.Matrix=makehgtform('translate',[m_estim(1 + 2*i) m_estim(1 + 2*i+1) 0],'zrotate',atan2(mi_eigenvec_sorted(2,1),mi_eigenvec_sorted(1,1)));
       
          line(gca,[ksi_groundtruth(1) m_groundtruth(1 + 2*i)],[ksi_groundtruth(2) m_groundtruth(1 + 2*i+1)],'LineStyle','--','Color',[0.75,0.75,0.75]);
       end
   hold off;
   
   refreshdata(h_fig);
   
   % save in each iteration
   dummy_eventdata.Key = 'space';
   keyInput(hObject,dummy_eventdata);
end

function keyInput(hObject,eventdata)
   global u_t;
   global z_mi;
   global h_timer;
   
   global m_groundtruth;
   global m_estim;
   global diffDriveModel;
   global ksi_estim;
   global ksi_groundtruth;
   global S_ksi_estim;
   global S_estim;
   global Q;
   global R;
   global time;

   %disp(eventdata.Key);
   
   switch eventdata.Key 
    case 'rightarrow';  u_t = u_t + [0; -0.1];
    case 'leftarrow'; u_t = u_t + [0; +0.1];
    case 'uparrow'; u_t = u_t + [+0.5; 0];
    case 'downarrow'; u_t = u_t + [-0.5; 0];
    
    case 'escape'; disp('Stopping h_timer / Deleting all timers'); stop(h_timer); listOfTimers = timerfindall; if ~isempty(listOfTimers); delete(listOfTimers(:)); end
    
    case 'space'
        load('SLAM_data');
        
        if (~exist('u_inputs','var'))
            u_inputs{1}.time = time;
            u_inputs{1}.u_t = u_t;
        else
            len = length(u_inputs);
            u_inputs{len+1}.time = time;
            u_inputs{len+1}.u_t = u_t;
        end
        
        if (~exist('m_estimates','var'))
            m_estimates{1}.time = time;
            m_estimates{1}.m_estim = m_estim;
        else
            len = length(m_estimates);
            m_estimates{len+1}.time = time;
            m_estimates{len+1}.m_estim = m_estim;
        end
       
        if (~exist('ksi_actual','var'))
            ksi_actual{1}.time = time;
            ksi_actual{1}.ksi_groundtruth = ksi_groundtruth;
        else
            len = length(ksi_actual);
            ksi_actual{len+1}.time = time;
            ksi_actual{len+1}.ksi_groundtruth = ksi_groundtruth;
        end
           
        if (~exist('ksi_estimates','var'))
            ksi_estimates{1}.time = time;
            ksi_estimates{1}.ksi_estim = ksi_estim;
        else
            len = length(ksi_estimates);
            ksi_estimates{len+1}.time = time;
            ksi_estimates{len+1}.ksi_estim = ksi_estim;
        end
        
        if (~exist('z_measurements','var'))
            z_measurements{1}.time = time;
            z_measurements{1}.z_mi = z_mi;
        else
            len = length(z_measurements);
            z_measurements{len+1}.time = time;
            z_measurements{len+1}.z_mi = z_mi;
        end
        
        save('SLAM_data','u_inputs','m_estimates','ksi_actual','ksi_estimates','z_measurements', '-append');
        
    otherwise  
   end
   
end

function closeFigureReq(src,callbackdata)
   % Close request function 
   % to display a question dialog box 
   %selection = questdlg('Close This Figure?',...
   %   'Close Request Function',...
   %   'Yes','No','Yes'); 
   %switch selection 
   %   case 'Yes'
   %      delete(gcf)
   %   case 'No'
   %   return 
   %end
   
   global h_timer;
   
   delete(src);

   disp('Stopping h_timer / Deleting all timers'); stop(h_timer); listOfTimers = timerfindall; if ~isempty(listOfTimers); delete(listOfTimers(:)); end
end
