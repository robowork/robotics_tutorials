clear all; clc;
global h_timer; %Handle for timer
global timer_period;
global time;

global val_f_x;
global val_f_y;
global val_f_z;
global val_m_x;
global val_m_y;
global val_m_z;
global scenario;
global fig_ax;
global h_fig_caz;
global h_fig_cel;
global sld_fx;
global sld_fy;
global sld_fz;
global sld_mx;
global sld_my;
global sld_mz;
global cbx_fx;
global cbx_fy;
global cbx_fz;
global cbx_mx;
global cbx_my;
global cbx_mz;
global rb_t0;
global ef_tx;
global ef_ty;
global ef_tz;
global ef_qw;
global ef_qx;
global ef_qy;
global ef_qz;
global S_p;
global S_R;
global S_T;
global B_p_x_0;
global B_p_y_0;
global B_p_z_0;
global B_p_qw_0;    
global B_p_qx_0;
global B_p_qy_0;
global B_p_qz_0; 
global omega_x_0;
global omega_y_0;
global omega_z_0;
global vel_x_0;
global vel_y_0;
global vel_z_0;

global B_p_x;
global B_p_y;
global B_p_z;
global B_p_qw;    
global B_p_qx;
global B_p_qy;
global B_p_qz; 
global omega_x;
global omega_y;
global omega_z;
global vel_x;
global vel_y;
global vel_z;
global dot_omega_x;
global dot_omega_y;
global dot_omega_z;
global dot_vel_x;
global dot_vel_y;
global dot_vel_z;

global ellipsoid_1_X;
global ellipsoid_1_Y;
global ellipsoid_1_Z;
global ellipsoid_1_axang;
global ellipsoid_2_X;
global ellipsoid_2_Y;
global ellipsoid_2_Z;
global ellipsoid_2_axang;
global combined_mass;
global combined_I;
bdclose all; close all;

timer_period = 0.20;
time = 0.0;

val_f_x = 0;
val_f_y = 0;
val_f_z = 0;
val_m_x = 0;
val_m_y = 0;
val_m_z = 0;

B_p_x_0 = 0.50;
B_p_y_0 = 0.25;
B_p_z_0 = 0.75;
B_p_qw_0 = 0.872;  
B_p_qx_0 = 0.215;
B_p_qy_0 = 0.189;
B_p_qz_0 = 0.398; 

omega_x_0 = 0*30;  %deg/s
omega_y_0 = 0*15;  %deg/s
omega_z_0 = 0*90;  %deg/s
vel_x_0 = 0*-0.5;
vel_y_0 = 0*-2.0;
vel_z_0 = 0*-0.25;

B_p_x = B_p_x_0;
B_p_y = B_p_y_0;
B_p_z = B_p_z_0;
B_p_qw = B_p_qw_0;  
B_p_qx = B_p_qx_0;
B_p_qy = B_p_qy_0;
B_p_qz = B_p_qz_0; 
omega_x = omega_x_0;  
omega_y = omega_y_0;
omega_z = omega_z_0;
vel_x = vel_x_0;
vel_y = vel_y_0;
vel_z = vel_z_0;
dot_omega_x = 0;  
dot_omega_y = 0;
dot_omega_z = 0;
dot_vel_x = 0;
dot_vel_y = 0;
dot_vel_z = 0;

scenario = 0;

% rigid body component (emulated as ellipsoids) masses and semiaxes
ellipsoid_1_mass = 1.0;
ellipsoid_1_a = 0.25;
ellipsoid_1_b = 0.75;
ellipsoid_1_c = 0.15;
ellipsoid_1_angle_z = deg2rad(45);
ellipsoid_1_angle_y = deg2rad(30);
ellipsoid_1_angle_x = deg2rad(15);

ellipsoid_2_mass = 5.0;
ellipsoid_2_a = 1.0;
ellipsoid_2_b = 0.10;
ellipsoid_2_c = 0.50;
ellipsoid_2_angle_z = deg2rad(-30);
ellipsoid_2_angle_y = deg2rad(-15);
ellipsoid_2_angle_x = deg2rad(45);

[ellipsoid_1_X,ellipsoid_1_Y,ellipsoid_1_Z] = ellipsoid(0,0,0,ellipsoid_1_a,ellipsoid_1_b,ellipsoid_1_c,10);
ellipsoid_1_I = (ellipsoid_1_mass/5)*[ellipsoid_1_b^2+ellipsoid_1_c^2 0                               0;
                                      0                               ellipsoid_1_a^2+ellipsoid_1_c^2 0;
                                      0                               0                               ellipsoid_1_a^2+ellipsoid_1_b^2];
ellipsoid_1_R_z = [cos(ellipsoid_1_angle_z) -sin(ellipsoid_1_angle_z) 0;
                    sin(ellipsoid_1_angle_z) cos(ellipsoid_1_angle_z)  0;
                    0            0             1];
ellipsoid_1_R_y = [cos(ellipsoid_1_angle_y)  0 sin(ellipsoid_1_angle_y);
                   0             1 0;
                   -sin(ellipsoid_1_angle_y) 0 cos(ellipsoid_1_angle_y)];
ellipsoid_1_R_x = [1 0            0;
                   0 cos(ellipsoid_1_angle_x) -sin(ellipsoid_1_angle_x);
                   0 sin(ellipsoid_1_angle_x) cos(ellipsoid_1_angle_x)];
ellipsoid_1_R_zyx = (ellipsoid_1_R_z * ellipsoid_1_R_y) * ellipsoid_1_R_x;
ellipsoid_1_axang = rotm2axang(ellipsoid_1_R_zyx);
ellipsoid_1_I = ellipsoid_1_R_zyx * ellipsoid_1_I * transpose(ellipsoid_1_R_zyx);

[ellipsoid_2_X,ellipsoid_2_Y,ellipsoid_2_Z] = ellipsoid(0,0,0,ellipsoid_2_a,ellipsoid_2_b,ellipsoid_2_c,10);
ellipsoid_2_I = (ellipsoid_2_mass/5)*[ellipsoid_2_b^2+ellipsoid_2_c^2 0                               0;
                                      0                               ellipsoid_2_a^2+ellipsoid_2_c^2 0;
                                      0                               0                               ellipsoid_2_a^2+ellipsoid_2_b^2];
ellipsoid_2_R_z = [cos(ellipsoid_2_angle_z) -sin(ellipsoid_2_angle_z) 0;
                    sin(ellipsoid_2_angle_z) cos(ellipsoid_2_angle_z)  0;
                    0            0             1];
ellipsoid_2_R_y = [cos(ellipsoid_2_angle_y)  0 sin(ellipsoid_2_angle_y);
                   0             1 0;
                   -sin(ellipsoid_2_angle_y) 0 cos(ellipsoid_2_angle_y)];
ellipsoid_2_R_x = [1 0            0;
                   0 cos(ellipsoid_2_angle_x) -sin(ellipsoid_2_angle_x);
                   0 sin(ellipsoid_2_angle_x) cos(ellipsoid_2_angle_x)];
ellipsoid_2_R_zyx = (ellipsoid_2_R_z * ellipsoid_2_R_y) * ellipsoid_2_R_x;
ellipsoid_2_axang = rotm2axang(ellipsoid_2_R_zyx);
ellipsoid_2_I = ellipsoid_2_R_zyx * ellipsoid_2_I * transpose(ellipsoid_2_R_zyx);

% combined rigid body mass and inertia tensor
combined_mass = ellipsoid_1_mass + ellipsoid_2_mass;
combined_I = ellipsoid_1_I + ellipsoid_2_I;

h_uifig = uifigure('Position',[10 10 550 850], 'Color','w', 'CloseRequestFcn',@closeFigureReq); %,'ValueChangedFcn',@updateFigure);
sld_fx = uislider('Parent',h_uifig,'Position',[100 825-0*50+0*5 400 10],'Limits',[-10 10],'Value',val_f_x);
sld_fy = uislider('Parent',h_uifig,'Position',[100 825-1*50+1*5 400 10],'Limits',[-10 10],'Value',val_f_y);
sld_fz = uislider('Parent',h_uifig,'Position',[100 825-2*50+2*5 400 10],'Limits',[-10 10],'Value',val_f_z);
sld_mx = uislider('Parent',h_uifig,'Position',[100 825-3*50+3*5 400 10],'Limits',[-1.0 1.0],'Value',val_m_x, 'Enable',1);
sld_my = uislider('Parent',h_uifig,'Position',[100 825-4*50+4*5 400 10],'Limits',[-1.0 1.0],'Value',val_m_y, 'Enable',1);
sld_mz = uislider('Parent',h_uifig,'Position',[100 825-5*50+5*5 400 10],'Limits',[-1.0 1.0],'Value',val_m_z, 'Enable',1);
fig_ax = uiaxes('Parent',h_uifig,'Position',[10 10 500 500]);
sld_fx.ValueChangingFcn = @updateUiSliderFx;
sld_fy.ValueChangingFcn = @updateUiSliderFy;
sld_fz.ValueChangingFcn = @updateUiSliderFz;
sld_mx.ValueChangingFcn = @updateUiSliderMx;
sld_my.ValueChangingFcn = @updateUiSliderMy;
sld_mz.ValueChangingFcn = @updateUiSliderMz;
cbx_fx = uicheckbox('Parent',h_uifig,'Position',[5 825-0*50+0*5 400 10],'Text','reset       fx','Value',0,...
                    'ValueChangedFcn',@(cbx_fx,event) updateUiCheckboxFx(cbx_fx,sld_fx));
cbx_fy = uicheckbox('Parent',h_uifig,'Position',[5 825-1*50+1*5 400 10],'Text','reset       fy','Value',0,...
                    'ValueChangedFcn',@(cbx_fy,event) updateUiCheckboxFy(cbx_fy,sld_fy));
cbx_fz = uicheckbox('Parent',h_uifig,'Position',[5 825-2*50+2*5 400 10],'Text','reset       fz','Value',0,...
                    'ValueChangedFcn',@(cbx_fz,event) updateUiCheckboxFz(cbx_fz,sld_fz));
cbx_mx = uicheckbox('Parent',h_uifig,'Position',[5 825-3*50+3*5 400 10],'Text','reset       mx','Value',0, 'Enable',1,...
                    'ValueChangedFcn',@(cbx_mx,event) updateUiCheckboxMx(cbx_mx,sld_mx));
cbx_my = uicheckbox('Parent',h_uifig,'Position',[5 825-4*50+4*5 400 10],'Text','reset       my','Value',0, 'Enable',1,...
                    'ValueChangedFcn',@(cbx_my,event) updateUiCheckboxMy(cbx_my,sld_my));
cbx_mz = uicheckbox('Parent',h_uifig,'Position',[5 825-5*50+5*5 400 10],'Text','reset       mz','Value',0, 'Enable',1,...
                    'ValueChangedFcn',@(cbx_mz,event) updateUiCheckboxMz(cbx_mz,sld_mz));
bg = uibuttongroup(h_uifig,'Position',[5 825-6*50+0*5 550-10 40], 'Backgroundcolor','w');
bg.SelectionChangedFcn = @updateUiRadiobutton;
rb_t0 = uiradiobutton('Parent',bg,'Position',[5+0*75 0 100 40],'Text','Wrench #1');
txa_t = uitextarea('Parent',h_uifig,'Position',[160 825-6*50+0*5+10 50 20],'Value','t_x,y,z','Editable','off','BackgroundColor',[0.75 0.75 0.75]);
ef_tx = uieditfield(h_uifig,'numeric','Position',[160+5+1*40 825-6*50+0*5+10 27.5 20],'Value',B_p_x_0,'ValueChangedFcn',@(ef_tx,event) updateUiTextboxTx(ef_tx,event));
ef_ty = uieditfield(h_uifig,'numeric','Position',[160+5+2*40 825-6*50+0*5+10 27.5 20],'Value',B_p_y_0,'ValueChangedFcn',@(ef_ty,event) updateUiTextboxTy(ef_ty,event));
ef_tz = uieditfield(h_uifig,'numeric','Position',[160+5+3*40 825-6*50+0*5+10 27.5 20],'Value',B_p_z_0,'ValueChangedFcn',@(ef_tz,event) updateUiTextboxTz(ef_tz,event));
txa_q = uitextarea('Parent',h_uifig,'Position',[160+5+4*40 825-6*50+0*5+10 75 20],'Value','q_w,x,y,z','Editable','off','BackgroundColor',[0.75 0.75 0.75]);
ef_qw = uieditfield(h_uifig,'numeric','Position',[160+20+5*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qw_0,'ValueChangedFcn',@(ef_qw,event) updateUiTextboxQw(ef_qw,event));
ef_qx = uieditfield(h_uifig,'numeric','Position',[160+20+6*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qx_0,'ValueChangedFcn',@(ef_qx,event) updateUiTextboxQx(ef_qx,event));
ef_qy = uieditfield(h_uifig,'numeric','Position',[160+20+7*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qy_0,'ValueChangedFcn',@(ef_qy,event) updateUiTextboxQy(ef_qy,event));
ef_qz = uieditfield(h_uifig,'numeric','Position',[160+20+8*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qz_0,'ValueChangedFcn',@(ef_qz,event) updateUiTextboxQz(ef_qz,event));

h_timer = timer('ExecutionMode','fixedRate','Period',timer_period,'TimerFcn',@updateTimer);

%%% Space Frame / origin %%%
S_p = [0;
       0;
       0];
S_R = [1 0 0;
       0 1 0;
       0 0 1];
S_T = [S_R    S_p;
       0 0 0  1];
   
translations = S_p';
rotations = rotm2quat(S_R);

plotTransforms(translations,rotations,'FrameSize',1.0,'Parent',fig_ax);
set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-3 3],'ylim',[-3 3],'zlim',[-3 3]);
dummy_hObject.SelectedObject = rb_t0;
updateUiRadiobutton(dummy_hObject,NaN);

[h_fig_caz,h_fig_cel] = view(fig_ax);

start(h_timer);

function updateTimer(hObject, eventdata)
   global h_timer; %Handle for timer
   global timer_period;
   global time;
   
   global B_p_x_0;
   global B_p_y_0;
   global B_p_z_0;
   global B_p_qw_0;    
   global B_p_qx_0;
   global B_p_qy_0;
   global B_p_qz_0;
   global omega_x_0;
   global omega_y_0;
   global omega_z_0;
   global vel_x_0;
   global vel_y_0;
   global vel_z_0;
   
   global B_p_x;
   global B_p_y;
   global B_p_z;
   global B_p_qw;    
   global B_p_qx;
   global B_p_qy;
   global B_p_qz;
   global omega_x;
   global omega_y;
   global omega_z;
   global vel_x;
   global vel_y;
   global vel_z;
   global dot_omega_x;
   global dot_omega_y;
   global dot_omega_z;
   global dot_vel_x;
   global dot_vel_y;
   global dot_vel_z;
   
   time = time + timer_period;
   if time > 5
       time = 0;
       
       B_p_x = B_p_x_0;
       B_p_y = B_p_y_0;
       B_p_z = B_p_z_0;
       B_p_qw = B_p_qw_0;  
       B_p_qx = B_p_qx_0;
       B_p_qy = B_p_qy_0;
       B_p_qz = B_p_qz_0; 
       omega_x = omega_x_0;  
       omega_y = omega_y_0;
       omega_z = omega_z_0;
       vel_x = vel_x_0;
       vel_y = vel_y_0;
       vel_z = vel_z_0;
       dot_omega_x = 0;  
       dot_omega_y = 0;
       dot_omega_z = 0;
       dot_vel_x = 0;
       dot_vel_y = 0;
       dot_vel_z = 0;
   end
   
   updateFigure();
end

function updateFigure()
   global h_timer; %Handle for timer
   global timer_period;
   global time;
  
   global val_f_x;
   global val_f_y;
   global val_f_z;
   global val_m_x;
   global val_m_y;
   global val_m_z;
   global scenario;
   global fig_ax;
   global h_fig_caz;
   global h_fig_cel;
   global S_p;
   global S_R;
   global S_T;
   
   global B_p_x;
   global B_p_y;
   global B_p_z;
   global B_p_qw;    
   global B_p_qx;
   global B_p_qy;
   global B_p_qz;
   global omega_x;
   global omega_y;
   global omega_z;
   global vel_x;
   global vel_y;
   global vel_z;
   global dot_omega_x;
   global dot_omega_y;
   global dot_omega_z;
   global dot_vel_x;
   global dot_vel_y;
   global dot_vel_z;

   global ellipsoid_1_X;
   global ellipsoid_1_Y;
   global ellipsoid_1_Z;
   global ellipsoid_1_axang;
   global ellipsoid_2_X;
   global ellipsoid_2_Y;
   global ellipsoid_2_Z;
   global ellipsoid_2_axang;
   global combined_mass;
   global combined_I;

   [h_fig_caz,h_fig_cel] = view(fig_ax);
   
   %%% Body Frame %%% 

   switch scenario 
    case 0
      % Twist expressed w.r.t. Body Frame
      B_p_zyx = [B_p_x;
                 B_p_y;
                 B_p_z];
      B_T_zyx_trans = [1 0 0  B_p_zyx(1);
                       0 1 0  B_p_zyx(2);
                       0 0 1  B_p_zyx(3);
                       0 0 0  1];
%       angle_z = deg2rad(45);
%       angle_y = deg2rad(30);
%       angle_x = deg2rad(15);
%       B_R_z = [cos(angle_z) -sin(angle_z) 0;
%                sin(angle_z) cos(angle_z)  0;
%                0            0             1];
%       B_T_z = [B_R_z  [0;0;0];
%                0 0 0  1];
%       B_R_y = [cos(angle_y)  0 sin(angle_y);
%                0             1 0;
%                -sin(angle_y) 0 cos(angle_y)];
%       B_T_y = [B_R_y  [0;0;0];
%                0 0 0  1];
%       B_R_x = [1 0            0;
%                0 cos(angle_x) -sin(angle_x);
%                0 sin(angle_x) cos(angle_x)];
%       B_T_x = [B_R_x  [0;0;0];
%                0 0 0  1];
%       B_T_zyx = (((S_T  *  B_T_zyx_trans) * B_T_z) * B_T_y) * B_T_x;
      B_T_zyx_rot = [quat2rotm([B_p_qw B_p_qx B_p_qy B_p_qz]) zeros(3,1);
                     zeros(1,3)                               1];
      B_T_zyx = (S_T  *  B_T_zyx_trans) * B_T_zyx_rot;          
      
      %%% Current Twist - Use to update KINEMATIC state %%%
      Twist_B = [deg2rad(omega_x);
                 deg2rad(omega_y);
                 deg2rad(omega_z);
                 vel_x;
                 vel_y;
                 vel_z];
        
      if norm(Twist_B(1:3)) < 1e-6
        V = eye(3);
        e_VMatrix_theta = [eye(3)     V*Twist_B(4:6) * timer_period;
                           zeros(1,3) 1];
      else
        omegaSkew = [0           -Twist_B(3) Twist_B(2);
                     Twist_B(3)  0           -Twist_B(1);
                     -Twist_B(2) Twist_B(1)  0] * timer_period;
        theta = norm(Twist_B(1:3) * timer_period);
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_VMatrix_theta = [e_omegaSkew V*Twist_B(4:6) * timer_period;
                           zeros(1,3)  1];
      end

      % Forward kinematics - compute updated pose
      B_T = B_T_zyx * e_VMatrix_theta;  % Post-multiply
      
      B_p_x = B_T(1,4);
      B_p_y = B_T(2,4);
      B_p_z = B_T(3,4);
      B_p_q = rotm2quat(B_T(1:3,1:3));
      B_p_qw = B_p_q(1);
      B_p_qx = B_p_q(2);
      B_p_qy = B_p_q(3);
      B_p_qz = B_p_q(4);
                   
      %%% Current Wrench - Use to update DYNAMIC state %%%
      Wrench_B = [val_m_x;
                  val_m_y;
                  val_m_z;
                  val_f_x;
                  val_f_y;
                  val_f_z];

      % 6x6 Inertia
      G = [combined_I zeros(3,3);
           zeros(3,3) combined_mass* eye(3)];
      % Twist
      omega_B_Skew = [0           -Twist_B(3) Twist_B(2);
                      Twist_B(3)  0           -Twist_B(1);
                      -Twist_B(2) Twist_B(1)  0];
      v_B_Skew = [0           -Twist_B(6) Twist_B(5);
                  Twist_B(6)  0           -Twist_B(4);
                  -Twist_B(5) Twist_B(4)  0];
      adv_B_transp = [omega_B_Skew v_B_Skew;
                      zeros(3,3)   omega_B_Skew];
      
      % Forward dynamics - compute updated twist
      dot_Twist_B = G \ (Wrench_B + adv_B_transp * G * Twist_B);
      
      Twist_B = Twist_B + dot_Twist_B * timer_period;
      omega_x = rad2deg(Twist_B(1));
      omega_y = rad2deg(Twist_B(2));
      omega_z = rad2deg(Twist_B(3));
      vel_x = Twist_B(4);
      vel_y = Twist_B(5);
      vel_z = Twist_B(6);

      %%% Visualize %%%
      translations = [S_T(1:3,4)';
                      B_T_zyx(1:3,4)';
                      B_T(1:3,4)'];
      rotations = [rotm2quat(S_R(1:3,1:3));
                   rotm2quat(B_T_zyx(1:3,1:3));
                   rotm2quat(B_T(1:3,1:3))];
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.25,'Parent',fig_ax);
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',0.5,'Parent',fig_ax);      
      plotTransforms(translations(3,:),rotations(3,:),'FrameSize',1.0,'Parent',fig_ax); 
      
      combined_trans = B_T(1:3,4);
      combined_axang = rotm2axang(B_T(1:3,1:3));
      hold(fig_ax,'on');
      h_ellipsoid_1 = surf(fig_ax,ellipsoid_1_X,ellipsoid_1_Y,ellipsoid_1_Z);
      rotate(h_ellipsoid_1,ellipsoid_1_axang(1:3),rad2deg(ellipsoid_1_axang(4)));
      h_ellipsoid_2 = surf(fig_ax,ellipsoid_2_X,ellipsoid_2_Y,ellipsoid_2_Z);
      rotate(h_ellipsoid_2,ellipsoid_2_axang(1:3),rad2deg(ellipsoid_2_axang(4)));
      hold(fig_ax,'off');
      rotate(h_ellipsoid_1,combined_axang(1:3),rad2deg(combined_axang(4)));
      rotate(h_ellipsoid_2,combined_axang(1:3),rad2deg(combined_axang(4)));
      h_ellipsoid_1_xdata = get(h_ellipsoid_1,'XData');  set(h_ellipsoid_1,'XData',h_ellipsoid_1_xdata+combined_trans(1));
      h_ellipsoid_1_ydata = get(h_ellipsoid_1,'YData');  set(h_ellipsoid_1,'YData',h_ellipsoid_1_ydata+combined_trans(2));
      h_ellipsoid_1_zdata = get(h_ellipsoid_1,'ZData');  set(h_ellipsoid_1,'ZData',h_ellipsoid_1_zdata+combined_trans(3));
      h_ellipsoid_2_xdata = get(h_ellipsoid_2,'XData');  set(h_ellipsoid_2,'XData',h_ellipsoid_2_xdata+combined_trans(1));
      h_ellipsoid_2_ydata = get(h_ellipsoid_2,'YData');  set(h_ellipsoid_2,'YData',h_ellipsoid_2_ydata+combined_trans(2));
      h_ellipsoid_2_zdata = get(h_ellipsoid_2,'ZData');  set(h_ellipsoid_2,'ZData',h_ellipsoid_2_zdata+combined_trans(3));
      
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-3 3],'ylim',[-3 3],'zlim',[-3 3]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Rigid Body dynamics given frame and motion initial conditions (Body Wrench)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    otherwise  
      bdclose all; close all;
 
   end
end

function updateUiCheckboxFx(hObject,sldObject)
   global val_f_x;
   
   if hObject.Value
       val_f_x = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxFy(hObject,sldObject)
   global val_f_y;
   
   if hObject.Value
       val_f_y = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxFz(hObject,sldObject)
   global val_f_z;
   
   if hObject.Value
       val_f_z = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxMx(hObject,sldObject)
   global val_m_x;
   
   if hObject.Value
       val_m_x = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxMy(hObject,sldObject)
   global val_m_y;
   
   if hObject.Value
       val_m_y = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxMz(hObject,sldObject)
   global val_m_z;
   
   if hObject.Value
       val_m_z = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end

function updateUiSliderFx(hObject,eventdata)
   global val_f_x;

   val_f_x = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderFy(hObject,eventdata)
   global val_f_y;

   val_f_y = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderFz(hObject,eventdata)
   global val_f_z;

   val_f_z = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderMx(hObject,eventdata)
   global val_m_x;

   val_m_x = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderMy(hObject,eventdata)
   global val_m_y;

   val_m_y = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderMz(hObject,eventdata)
   global val_m_z;

   val_m_z = eventdata.Value;
   
   %updateFigure();
end

function updateUiTextboxTx(hObject,eventdata)
   global B_p_x_0;
   
   B_p_x_0 = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxTy(hObject,eventdata)
   global B_p_y_0;
   
   B_p_y_0 = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxTz(hObject,eventdata)
   global B_p_z_0;
   
   B_p_z_0 = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQw(hObject,eventdata)
   global B_p_qw_0;
   
   B_p_qw_0 = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQx(hObject,eventdata)
   global B_p_qx_0;
   
   B_p_qx_0 = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQy(hObject,eventdata)
   global B_p_qy_0;
   
   B_p_qy_0 = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQz(hObject,eventdata)
   global B_p_qz_0;
   
   B_p_qz_0 = eventdata.Value;

   %updateFigure();
end

function updateUiRadiobutton(hObject,eventdata)
   global val_f_x;
   global val_f_y;
   global val_f_z;
   global val_m_x;
   global val_m_y;
   global val_m_z;
   global scenario;
   global sld_fx;
   global sld_fy;
   global sld_fz;
   global sld_mx;
   global sld_my;
   global sld_mz;
   global cbx_fx;
   global cbx_fy;
   global cbx_fz;
   global cbx_mx;
   global cbx_my;
   global cbx_mz;   
   global rb_t0;
   
   val_f_x = 0; sld_fx.Value = val_f_x; 
   val_f_y = 0; sld_fy.Value = val_f_y; 
   val_f_z = 0; sld_fz.Value = val_f_z; 
   val_m_x = 0; sld_mx.Value = val_m_x; 
   val_m_y = 0; sld_my.Value = val_m_y; 
   val_m_z = 0; sld_mz.Value = val_m_z; 
           
   switch hObject.SelectedObject 
     case rb_t0
       scenario =  0;
   end
   
   updateFigure();
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