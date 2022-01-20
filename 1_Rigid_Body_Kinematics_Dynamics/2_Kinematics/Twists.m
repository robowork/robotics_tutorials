clear all; clc;
global h_timer; %Handle for timer
global timer_period;
global time;

global val_omega_x;
global val_omega_y;
global val_omega_z;
global val_vel_x;
global val_vel_y;
global val_vel_z;
global scenario;
global fig_ax;
global h_fig_caz;
global h_fig_cel;
global sld_wx;
global sld_wy;
global sld_wz;
global sld_vx;
global sld_vy;
global sld_vz;
global cbx_wx;
global cbx_wy;
global cbx_wz;
global cbx_vx;
global cbx_vy;
global cbx_vz;
global rb_t0;
global rb_t1;
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
global B_p_x;
global B_p_y;
global B_p_z;
global B_p_qw;    
global B_p_qx;
global B_p_qy;
global B_p_qz; 
bdclose all; close all;

timer_period = 0.20;
time = 0.0;

val_omega_x = 0;
val_omega_y = 0;
val_omega_z = 0;
val_vel_x = 0;
val_vel_y = 0;
val_vel_z = 0;

B_p_x = 0.50;
B_p_y = 0.25;
B_p_z = 0.75;
B_p_qw = 0.872;  
B_p_qx = 0.215;
B_p_qy = 0.189;
B_p_qz = 0.398; 

scenario = 0;

h_uifig = uifigure('Position',[10 10 550 850], 'Color','w', 'CloseRequestFcn',@closeFigureReq); %,'ValueChangedFcn',@updateFigure);
sld_wx = uislider('Parent',h_uifig,'Position',[100 825-0*50+0*5 400 10],'Limits',[-180 180],'Value',val_omega_x);
sld_wy = uislider('Parent',h_uifig,'Position',[100 825-1*50+1*5 400 10],'Limits',[-180 180],'Value',val_omega_y);
sld_wz = uislider('Parent',h_uifig,'Position',[100 825-2*50+2*5 400 10],'Limits',[-180 180],'Value',val_omega_z);
sld_vx = uislider('Parent',h_uifig,'Position',[100 825-3*50+3*5 400 10],'Limits',[-1.0 1.0],'Value',val_vel_x, 'Enable',1);
sld_vy = uislider('Parent',h_uifig,'Position',[100 825-4*50+4*5 400 10],'Limits',[-1.0 1.0],'Value',val_vel_y, 'Enable',1);
sld_vz = uislider('Parent',h_uifig,'Position',[100 825-5*50+5*5 400 10],'Limits',[-1.0 1.0],'Value',val_vel_z, 'Enable',1);
fig_ax = uiaxes('Parent',h_uifig,'Position',[10 10 500 500]);
sld_wx.ValueChangingFcn = @updateUiSliderWx;
sld_wy.ValueChangingFcn = @updateUiSliderWy;
sld_wz.ValueChangingFcn = @updateUiSliderWz;
sld_vx.ValueChangingFcn = @updateUiSliderVx;
sld_vy.ValueChangingFcn = @updateUiSliderVy;
sld_vz.ValueChangingFcn = @updateUiSliderVz;
cbx_wx = uicheckbox('Parent',h_uifig,'Position',[5 825-0*50+0*5 400 10],'Text','reset       ωx','Value',0,...
                    'ValueChangedFcn',@(cbx_wx,event) updateUiCheckboxWx(cbx_wx,sld_wx));
cbx_wy = uicheckbox('Parent',h_uifig,'Position',[5 825-1*50+1*5 400 10],'Text','reset       ωy','Value',0,...
                    'ValueChangedFcn',@(cbx_wy,event) updateUiCheckboxWy(cbx_wy,sld_wy));
cbx_wz = uicheckbox('Parent',h_uifig,'Position',[5 825-2*50+2*5 400 10],'Text','reset       ωz','Value',0,...
                    'ValueChangedFcn',@(cbx_wz,event) updateUiCheckboxWz(cbx_wz,sld_wz));
cbx_vx = uicheckbox('Parent',h_uifig,'Position',[5 825-3*50+3*5 400 10],'Text','reset       vx','Value',0, 'Enable',1,...
                    'ValueChangedFcn',@(cbx_vx,event) updateUiCheckboxVx(cbx_vx,sld_vx));
cbx_vy = uicheckbox('Parent',h_uifig,'Position',[5 825-4*50+4*5 400 10],'Text','reset       vy','Value',0, 'Enable',1,...
                    'ValueChangedFcn',@(cbx_vy,event) updateUiCheckboxVy(cbx_vy,sld_vy));
cbx_vz = uicheckbox('Parent',h_uifig,'Position',[5 825-5*50+5*5 400 10],'Text','reset       vz','Value',0, 'Enable',1,...
                    'ValueChangedFcn',@(cbx_vz,event) updateUiCheckboxVz(cbx_vz,sld_vz));
bg = uibuttongroup(h_uifig,'Position',[5 825-6*50+0*5 550-10 40], 'Backgroundcolor','w');
bg.SelectionChangedFcn = @updateUiRadiobutton;
rb_t0 = uiradiobutton('Parent',bg,'Position',[5+0*75 0 100 40],'Text','Twist #1');
rb_t1 = uiradiobutton('Parent',bg,'Position',[5+1*75 0 100 40],'Text','Twist #2');
txa_t = uitextarea('Parent',h_uifig,'Position',[160 825-6*50+0*5+10 50 20],'Value','t_x,y,z','Editable','off','BackgroundColor',[0.75 0.75 0.75]);
ef_tx = uieditfield(h_uifig,'numeric','Position',[160+5+1*40 825-6*50+0*5+10 27.5 20],'Value',B_p_x,'ValueChangedFcn',@(ef_tx,event) updateUiTextboxTx(ef_tx,event));
ef_ty = uieditfield(h_uifig,'numeric','Position',[160+5+2*40 825-6*50+0*5+10 27.5 20],'Value',B_p_y,'ValueChangedFcn',@(ef_ty,event) updateUiTextboxTy(ef_ty,event));
ef_tz = uieditfield(h_uifig,'numeric','Position',[160+5+3*40 825-6*50+0*5+10 27.5 20],'Value',B_p_z,'ValueChangedFcn',@(ef_tz,event) updateUiTextboxTz(ef_tz,event));
txa_q = uitextarea('Parent',h_uifig,'Position',[160+5+4*40 825-6*50+0*5+10 75 20],'Value','q_w,x,y,z','Editable','off','BackgroundColor',[0.75 0.75 0.75]);
ef_qw = uieditfield(h_uifig,'numeric','Position',[160+20+5*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qw,'ValueChangedFcn',@(ef_qw,event) updateUiTextboxQw(ef_qw,event));
ef_qx = uieditfield(h_uifig,'numeric','Position',[160+20+6*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qx,'ValueChangedFcn',@(ef_qx,event) updateUiTextboxQx(ef_qx,event));
ef_qy = uieditfield(h_uifig,'numeric','Position',[160+20+7*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qy,'ValueChangedFcn',@(ef_qy,event) updateUiTextboxQy(ef_qy,event));
ef_qz = uieditfield(h_uifig,'numeric','Position',[160+20+8*40 825-6*50+0*5+10 37.5 20],'Value',B_p_qz,'ValueChangedFcn',@(ef_qz,event) updateUiTextboxQz(ef_qz,event));

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
   
   time = time + timer_period;
   if time > 2
       time = 0;
   end
   
   updateFigure();
end

function updateFigure()
   global h_timer; %Handle for timer
   global timer_period;
   global time;
  
   global val_omega_x;
   global val_omega_y;
   global val_omega_z;
   global val_vel_x;
   global val_vel_y;
   global val_vel_z;
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
      
      Twist_B = [deg2rad(val_omega_x);
                 deg2rad(val_omega_y);
                 deg2rad(val_omega_z);
                 val_vel_x;
                 val_vel_y;
                 val_vel_z];
        
      if time < 1e-3
        e_VMatrix_theta = [eye(3)     zeros(3,1);
                           zeros(1,3) 1];  
      elseif norm(Twist_B(1:3)) < 1e-6
        V = eye(3);
        e_VMatrix_theta = [eye(3)     V*Twist_B(4:6) * time;
                           zeros(1,3) 1];
      else
        omegaSkew = [0           -Twist_B(3) Twist_B(2);
                     Twist_B(3)  0           -Twist_B(1);
                     -Twist_B(2) Twist_B(1)  0] * time;
        theta = norm(Twist_B(1:3) * time);
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_VMatrix_theta = [e_omegaSkew V*Twist_B(4:6) * time;
                           zeros(1,3)  1];
      end

      B_T = B_T_zyx * e_VMatrix_theta;  % Post-multiply
      
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
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-3 3],'ylim',[-3 3],'zlim',[-3 3]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Screw Motion given frame initial conditions (Twist w.r.t. Body Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);

    case 1
      % Twist expressed w.r.t. Space Frame
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
      
      Twist_B = [deg2rad(val_omega_x);
                 deg2rad(val_omega_y);
                 deg2rad(val_omega_z);
                 val_vel_x;
                 val_vel_y;
                 val_vel_z];
        
      if time < 1e-3
        e_VMatrix_theta = [eye(3)     zeros(3,1);
                           zeros(1,3) 1];  
      elseif norm(Twist_B(1:3)) < 1e-6
        V = eye(3);
        e_VMatrix_theta = [eye(3)     V*Twist_B(4:6) * time;
                           zeros(1,3) 1];
      else
        omegaSkew = [0           -Twist_B(3) Twist_B(2);
                     Twist_B(3)  0           -Twist_B(1);
                     -Twist_B(2) Twist_B(1)  0] * time;
        theta = norm(Twist_B(1:3) * time);
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_VMatrix_theta = [e_omegaSkew V*Twist_B(4:6) * time;
                           zeros(1,3)  1];
      end

      B_T = e_VMatrix_theta * B_T_zyx;  % Pre-multiply
      
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
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-3 3],'ylim',[-3 3],'zlim',[-3 3]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Screw Motion given frame initial conditions (Twist w.r.t. Space Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    otherwise  
      bdclose all; close all;
 
   end
end

function updateUiCheckboxWx(hObject,sldObject)
   global val_omega_x;
   
   if hObject.Value
       val_omega_x = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxWy(hObject,sldObject)
   global val_omega_y;
   
   if hObject.Value
       val_omega_y = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxWz(hObject,sldObject)
   global val_omega_z;
   
   if hObject.Value
       val_omega_z = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxVx(hObject,sldObject)
   global val_vel_x;
   
   if hObject.Value
       val_vel_x = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxVy(hObject,sldObject)
   global val_vel_y;
   
   if hObject.Value
       val_vel_y = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxVz(hObject,sldObject)
   global val_vel_z;
   
   if hObject.Value
       val_vel_z = 0;

       %updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end

function updateUiSliderWx(hObject,eventdata)
   global val_omega_x;

   val_omega_x = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderWy(hObject,eventdata)
   global val_omega_y;

   val_omega_y = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderWz(hObject,eventdata)
   global val_omega_z;

   val_omega_z = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderVx(hObject,eventdata)
   global val_vel_x;

   val_vel_x = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderVy(hObject,eventdata)
   global val_vel_y;

   val_vel_y = eventdata.Value;
   
   %updateFigure();
end
function updateUiSliderVz(hObject,eventdata)
   global val_vel_z;

   val_vel_z = eventdata.Value;
   
   %updateFigure();
end

function updateUiTextboxTx(hObject,eventdata)
   global B_p_x;
   
   B_p_x = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxTy(hObject,eventdata)
   global B_p_y;
   
   B_p_y = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxTz(hObject,eventdata)
   global B_p_z;
   
   B_p_z = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQw(hObject,eventdata)
   global B_p_qw;
   
   B_p_qw = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQx(hObject,eventdata)
   global B_p_qx;
   
   B_p_qx = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQy(hObject,eventdata)
   global B_p_qy;
   
   B_p_qy = eventdata.Value;

   %updateFigure();
end
function updateUiTextboxQz(hObject,eventdata)
   global B_p_qz;
   
   B_p_qz = eventdata.Value;

   %updateFigure();
end

function updateUiRadiobutton(hObject,eventdata)
   global val_omega_x;
   global val_omega_y;
   global val_omega_z;
   global val_vel_x;
   global val_vel_y;
   global val_vel_z;
   global scenario;
   global sld_wx;
   global sld_wy;
   global sld_wz;
   global sld_vx;
   global sld_vy;
   global sld_vz;
   global cbx_wx;
   global cbx_wy;
   global cbx_wz;
   global cbx_vx;
   global cbx_vy;
   global cbx_vz;   
   global rb_t0;
   global rb_t1;
   
   val_omega_x = 0; sld_wx.Value = val_omega_x; 
   val_omega_y = 0; sld_wy.Value = val_omega_y; 
   val_omega_z = 0; sld_wz.Value = val_omega_z; 
   val_vel_x = 0; sld_vx.Value = val_vel_x; 
   val_vel_y = 0; sld_vy.Value = val_vel_y; 
   val_vel_z = 0; sld_vz.Value = val_vel_z; 
           
   switch hObject.SelectedObject 
     case rb_t0
       scenario =  0;
     case rb_t1
       scenario =  1;
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