clear all; clc;
global val_rot_x;
global val_rot_y;
global val_rot_z;
global val_trans_x;
global val_trans_y;
global val_trans_z;
global scenario;
global fig_ax;
global h_fig_caz;
global h_fig_cel;
global sld_rx;
global sld_ry;
global sld_rz;
global sld_px;
global sld_py;
global sld_pz;
global cbx_rx;
global cbx_ry;
global cbx_rz;
global cbx_px;
global cbx_py;
global cbx_pz;
global rb_r0;
global rb_r1;
global rb_r2;
global rb_p0;
global rb_p1;
global rb_p2;
global S_p;
global S_R;
global S_T;
bdclose all; close all;

val_rot_x = 0;
val_rot_y = 0;
val_rot_z = 0;
val_trans_x = 0;
val_trans_y = 0;
val_trans_z = 0;

scenario = 0;

h_uifig = uifigure('Position',[10 10 550 850], 'Color','w'); %,'ValueChangedFcn',@updateFigure);
sld_rx = uislider('Parent',h_uifig,'Position',[100 825-0*50+0*5 400 10],'Limits',[-180 180],'Value',val_rot_x);
sld_ry = uislider('Parent',h_uifig,'Position',[100 825-1*50+1*5 400 10],'Limits',[-180 180],'Value',val_rot_y);
sld_rz = uislider('Parent',h_uifig,'Position',[100 825-2*50+2*5 400 10],'Limits',[-180 180],'Value',val_rot_z);
sld_px = uislider('Parent',h_uifig,'Position',[100 825-3*50+3*5 400 10],'Limits',[-1 1],'Value',val_trans_x, 'Enable',0);
sld_py = uislider('Parent',h_uifig,'Position',[100 825-4*50+4*5 400 10],'Limits',[-1 1],'Value',val_trans_y, 'Enable',0);
sld_pz = uislider('Parent',h_uifig,'Position',[100 825-5*50+5*5 400 10],'Limits',[-1 1],'Value',val_trans_z, 'Enable',0);
fig_ax = uiaxes('Parent',h_uifig,'Position',[10 10 500 500]);
sld_rx.ValueChangingFcn = @updateUiSliderRx;
sld_ry.ValueChangingFcn = @updateUiSliderRy;
sld_rz.ValueChangingFcn = @updateUiSliderRz;
sld_px.ValueChangingFcn = @updateUiSliderPx;
sld_py.ValueChangingFcn = @updateUiSliderPy;
sld_pz.ValueChangingFcn = @updateUiSliderPz;
cbx_rx = uicheckbox('Parent',h_uifig,'Position',[5 825-0*50+0*5 400 10],'Text','reset       Rx','Value',0,...
                    'ValueChangedFcn',@(cbx_rx,event) updateUiCheckboxRx(cbx_rx,sld_rx));
cbx_ry = uicheckbox('Parent',h_uifig,'Position',[5 825-1*50+1*5 400 10],'Text','reset       Ry','Value',0,...
                    'ValueChangedFcn',@(cbx_ry,event) updateUiCheckboxRy(cbx_ry,sld_ry));
cbx_rz = uicheckbox('Parent',h_uifig,'Position',[5 825-2*50+2*5 400 10],'Text','reset       Rz','Value',0,...
                    'ValueChangedFcn',@(cbx_rz,event) updateUiCheckboxRz(cbx_rz,sld_rz));
cbx_px = uicheckbox('Parent',h_uifig,'Position',[5 825-3*50+3*5 400 10],'Text','reset       px','Value',0, 'Enable',0,...
                    'ValueChangedFcn',@(cbx_px,event) updateUiCheckboxPx(cbx_px,sld_px));
cbx_py = uicheckbox('Parent',h_uifig,'Position',[5 825-4*50+4*5 400 10],'Text','reset       py','Value',0, 'Enable',0,...
                    'ValueChangedFcn',@(cbx_py,event) updateUiCheckboxPy(cbx_py,sld_py));
cbx_pz = uicheckbox('Parent',h_uifig,'Position',[5 825-5*50+5*5 400 10],'Text','reset       pz','Value',0, 'Enable',0,...
                    'ValueChangedFcn',@(cbx_pz,event) updateUiCheckboxPz(cbx_pz,sld_pz));
bg = uibuttongroup(h_uifig,'Position',[5 825-6*50+0*5 550-10 40], 'Backgroundcolor','w');
bg.SelectionChangedFcn = @updateUiRadiobutton;
rb_r0 = uiradiobutton('Parent',bg,'Position',[5+0*75 0 100 40],'Text','Rot #1');
rb_r1 = uiradiobutton('Parent',bg,'Position',[5+1*75 0 100 40],'Text','Rot #2');
rb_r2 = uiradiobutton('Parent',bg,'Position',[5+2*75 0 100 40],'Text','Rot #3');
rb_p0 = uiradiobutton('Parent',bg,'Position',[5+3*75+1*10 0 100 40],'Text','Homog #1');
rb_p1 = uiradiobutton('Parent',bg,'Position',[5+4*75+2*10 0 100 40],'Text','Homog #2');
rb_p2 = uiradiobutton('Parent',bg,'Position',[5+5*75+3*10 0 100 40],'Text','Homog #3');

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
set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
dummy_hObject.SelectedObject = rb_r0;
updateUiRadiobutton(dummy_hObject,NaN);

[h_fig_caz,h_fig_cel] = view(fig_ax);

function updateFigure()
   global val_rot_x;
   global val_rot_y;
   global val_rot_z;
   global val_trans_x;
   global val_trans_y;
   global val_trans_z;
   global scenario;
   global fig_ax;
   global h_fig_caz;
   global h_fig_cel;
   global S_p;
   global S_R;
   global S_T;

   [h_fig_caz,h_fig_cel] = view(fig_ax);
   
   %%% Body Frame %%% 

   switch scenario 
    case 0
      % Rotation sequence around z-y-x-axes (post-multiply, each sequential rotation w.r.t. Body Frame)
      angle_z = deg2rad(val_rot_z);
      angle_y = deg2rad(val_rot_y);
      angle_x = deg2rad(val_rot_x);
      B_R_z = [cos(angle_z) -sin(angle_z) 0;
               sin(angle_z) cos(angle_z)  0;
               0            0             1];
      B_R_y = [cos(angle_y)  0 sin(angle_y);
               0             1 0;
               -sin(angle_y) 0 cos(angle_y)];
      B_R_x = [1 0            0;
               0 cos(angle_x) -sin(angle_x);
               0 sin(angle_x) cos(angle_x)];
        
      B_p_zyx = [0;
                 0;
                 0];
      B_R_zyx = ((S_R  *  B_R_z) * B_R_y) * B_R_x;

      translations = [S_p';
                      B_p_zyx'];
      rotations = [rotm2quat(S_R);
                   rotm2quat(B_R_zyx)];
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.25,'Parent',fig_ax);
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',1.0,'Parent',fig_ax);      
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Rotation sequence around z-y-x-axes';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    case 1
      % Post-multiply, rotation w.r.t. Body Frame
      angle_z = deg2rad(45);
      angle_y = deg2rad(30);
      angle_x = deg2rad(15);
      B_R_z = [cos(angle_z) -sin(angle_z) 0;
               sin(angle_z) cos(angle_z)  0;
               0            0             1];
      B_R_y = [cos(angle_y)  0 sin(angle_y);
               0             1 0;
               -sin(angle_y) 0 cos(angle_y)];
      B_R_x = [1 0            0;
               0 cos(angle_x) -sin(angle_x);
               0 sin(angle_x) cos(angle_x)];         
      B_p_zyx = [0;
                 0;
                 0];
      B_R_zyx = ((S_R  *  B_R_z) * B_R_y) * B_R_x;
      
      B_extra_R_z = axang2rotm([0 0 1 deg2rad(val_rot_z)]);
      B_extra_R_y = axang2rotm([0 1 0 deg2rad(val_rot_y)]);
      B_extra_R_x = axang2rotm([1 0 0 deg2rad(val_rot_x)]);
      
      B_p = B_p_zyx;
      B_R = ((B_R_zyx * B_extra_R_z) * B_extra_R_y) * B_extra_R_x;  % Post-multiply
      
      translations = [S_p';
                      B_p_zyx';
                      B_p'];
      rotations = [rotm2quat(S_R);
                   rotm2quat(B_R_zyx);
                   rotm2quat(B_R)];
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.25,'Parent',fig_ax);
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',0.5,'Parent',fig_ax);      
      plotTransforms(translations(3,:),rotations(3,:),'FrameSize',1.0,'Parent',fig_ax);      
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Post-multiply [ψ,θ,φ]=[45,30,15] (rotate w.r.t. Body Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);

    case 2
      % Pre-multiply, rotation w.r.t. Space Frame
      angle_z = deg2rad(45);
      angle_y = deg2rad(30);
      angle_x = deg2rad(15);
      B_R_z = [cos(angle_z) -sin(angle_z) 0;
               sin(angle_z) cos(angle_z)  0;
               0            0             1];
      B_R_y = [cos(angle_y)  0 sin(angle_y);
               0             1 0;
               -sin(angle_y) 0 cos(angle_y)];
      B_R_x = [1 0            0;
               0 cos(angle_x) -sin(angle_x);
               0 sin(angle_x) cos(angle_x)];
      B_p_zyx = [0;
                 0;
                 0];
      B_R_zyx = ((S_R  *  B_R_z) * B_R_y) * B_R_x;
      
      B_extra_R_z = axang2rotm([0 0 1 deg2rad(val_rot_z)]);
      B_extra_R_y = axang2rotm([0 1 0 deg2rad(val_rot_y)]);
      B_extra_R_x = axang2rotm([1 0 0 deg2rad(val_rot_x)]);
      
      B_p = B_p_zyx;
      B_R = B_extra_R_x * (B_extra_R_y * (B_extra_R_z * B_R_zyx));  % Pre-multiply
      
      translations = [S_p';
                      B_p_zyx';
                      B_p'];
      rotations = [rotm2quat(S_R);
                   rotm2quat(B_R_zyx);
                   rotm2quat(B_R)];
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.25,'Parent',fig_ax);
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',0.5,'Parent',fig_ax);      
      plotTransforms(translations(3,:),rotations(3,:),'FrameSize',1.0,'Parent',fig_ax);      
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Pre-multiply [ψ,θ,φ]=[45,30,15] (rotate w.r.t. Space Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    case 3
      % Translation w.r.t. all axes, Rotation sequence around z-y-x-axes (post-multiply, each sequential transformation w.r.t. Body Frame)
      B_p_zyx = [val_trans_x;
                 val_trans_y;
                 val_trans_z];
      B_T_zyx_trans = [1 0 0  B_p_zyx(1);
                       0 1 0  B_p_zyx(2);
                       0 0 1  B_p_zyx(3);
                       0 0 0  1];
      angle_z = deg2rad(val_rot_z);
      angle_y = deg2rad(val_rot_y);
      angle_x = deg2rad(val_rot_x);
      B_R_z = [cos(angle_z) -sin(angle_z) 0;
               sin(angle_z) cos(angle_z)  0;
               0            0             1];
      B_T_z = [B_R_z  [0;0;0];
               0 0 0  1];
      B_R_y = [cos(angle_y)  0 sin(angle_y);
               0             1 0;
               -sin(angle_y) 0 cos(angle_y)];
      B_T_y = [B_R_y  [0;0;0];
               0 0 0  1];
      B_R_x = [1 0            0;
               0 cos(angle_x) -sin(angle_x);
               0 sin(angle_x) cos(angle_x)];
      B_T_x = [B_R_x  [0;0;0];
               0 0 0  1];
           
      B_T_zyx = (((S_T  *  B_T_zyx_trans) * B_T_z) * B_T_y) * B_T_x;

      translations = [S_T(1:3,4)';
                      B_T_zyx(1:3,4)'];
      rotations = [rotm2quat(S_R(1:3,1:3));
                   rotm2quat(B_T_zyx(1:3,1:3))];
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.25,'Parent',fig_ax);
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',1.0,'Parent',fig_ax);      
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Translation w.r.t. all axes, Rotation sequence around z-y-x-axes';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    case 4
      % Post-multiply, transformation w.r.t. Body Frame
      B_p_zyx = [0.50;
                 0.25;
                 0.75];
      B_T_zyx_trans = [1 0 0  B_p_zyx(1);
                       0 1 0  B_p_zyx(2);
                       0 0 1  B_p_zyx(3);
                       0 0 0  1];
      angle_z = deg2rad(45);
      angle_y = deg2rad(30);
      angle_x = deg2rad(15);
      B_R_z = [cos(angle_z) -sin(angle_z) 0;
               sin(angle_z) cos(angle_z)  0;
               0            0             1];
      B_T_z = [B_R_z  [0;0;0];
               0 0 0  1];
      B_R_y = [cos(angle_y)  0 sin(angle_y);
               0             1 0;
               -sin(angle_y) 0 cos(angle_y)];
      B_T_y = [B_R_y  [0;0;0];
               0 0 0  1];
      B_R_x = [1 0            0;
               0 cos(angle_x) -sin(angle_x);
               0 sin(angle_x) cos(angle_x)];
      B_T_x = [B_R_x  [0;0;0];
               0 0 0  1];
      B_T_zyx = (((S_T  *  B_T_zyx_trans) * B_T_z) * B_T_y) * B_T_x;
      
      B_extra_p = [val_trans_x;
                   val_trans_y;
                   val_trans_z];
      B_extra_T_trans = [1 0 0  B_extra_p(1);
                         0 1 0  B_extra_p(2);
                         0 0 1  B_extra_p(3);
                         0 0 0  1];
      B_extra_R_z = axang2rotm([0 0 1 deg2rad(val_rot_z)]);
      B_extra_T_z = [B_extra_R_z  [0;0;0];
                     0 0 0  1];
      B_extra_R_y = axang2rotm([0 1 0 deg2rad(val_rot_y)]);
      B_extra_T_y = [B_extra_R_y  [0;0;0];
                     0 0 0  1];
      B_extra_R_x = axang2rotm([1 0 0 deg2rad(val_rot_x)]);
      B_extra_T_x = [B_extra_R_x  [0;0;0];
                     0 0 0  1];

      B_T = (((B_T_zyx * B_extra_T_trans) * B_extra_T_z) * B_extra_T_y) * B_extra_T_x;  % Post-multiply

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
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Post-multiply [x,y,z]=[0.5,0.25,0.75] [ψ,θ,φ]=[45,30,15] (transform w.r.t. Body Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);

    case 5
      % Pre-multiply, transformation w.r.t. Space Frame
      B_p_zyx = [0.50;
                 0.25;
                 0.75];
      B_T_zyx_trans = [1 0 0  B_p_zyx(1);
                       0 1 0  B_p_zyx(2);
                       0 0 1  B_p_zyx(3);
                       0 0 0  1];
      angle_z = deg2rad(45);
      angle_y = deg2rad(30);
      angle_x = deg2rad(15);
      B_R_z = [cos(angle_z) -sin(angle_z) 0;
               sin(angle_z) cos(angle_z)  0;
               0            0             1];
      B_T_z = [B_R_z  [0;0;0];
               0 0 0  1];
      B_R_y = [cos(angle_y)  0 sin(angle_y);
               0             1 0;
               -sin(angle_y) 0 cos(angle_y)];
      B_T_y = [B_R_y  [0;0;0];
               0 0 0  1];
      B_R_x = [1 0            0;
               0 cos(angle_x) -sin(angle_x);
               0 sin(angle_x) cos(angle_x)];
      B_T_x = [B_R_x  [0;0;0];
               0 0 0  1];
      B_T_zyx = (((S_T  *  B_T_zyx_trans) * B_T_z) * B_T_y) * B_T_x;
      
      B_extra_p = [val_trans_x;
                   val_trans_y;
                   val_trans_z];
      B_extra_T_trans = [1 0 0  B_extra_p(1);
                         0 1 0  B_extra_p(2);
                         0 0 1  B_extra_p(3);
                         0 0 0  1];
      B_extra_R_z = axang2rotm([0 0 1 deg2rad(val_rot_z)]);
      B_extra_T_z = [B_extra_R_z  [0;0;0];
                     0 0 0  1];
      B_extra_R_y = axang2rotm([0 1 0 deg2rad(val_rot_y)]);
      B_extra_T_y = [B_extra_R_y  [0;0;0];
                     0 0 0  1];
      B_extra_R_x = axang2rotm([1 0 0 deg2rad(val_rot_x)]);
      B_extra_T_x = [B_extra_R_x  [0;0;0];
                     0 0 0  1];

      B_T = B_extra_T_x * (B_extra_T_y * (B_extra_T_z * (B_extra_T_trans * B_T_zyx)));  % Pre-multiply

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
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-2 2],'ylim',[-2 2],'zlim',[-2 2]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = 'Pre-multiply [x,y,z]=[0.5,0.25,0.75] [ψ,θ,φ]=[45,30,15] (transform w.r.t. Space Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    otherwise  
      bdclose all; close all;
 
   end

end

function updateUiCheckboxRx(hObject,sldObject)
   global val_rot_x;
   
   if hObject.Value
       val_rot_x = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxRy(hObject,sldObject)
   global val_rot_y;
   
   if hObject.Value
       val_rot_y = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxRz(hObject,sldObject)
   global val_rot_z;
   
   if hObject.Value
       val_rot_z = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxPx(hObject,sldObject)
   global val_trans_x;
   
   if hObject.Value
       val_trans_x = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxPy(hObject,sldObject)
   global val_trans_y;
   
   if hObject.Value
       val_trans_y = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxPz(hObject,sldObject)
   global val_trans_z;
   
   if hObject.Value
       val_trans_z = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end

function updateUiSliderRx(hObject,eventdata)
   global val_rot_x;

   val_rot_x = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderRy(hObject,eventdata)
   global val_rot_y;

   val_rot_y = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderRz(hObject,eventdata)
   global val_rot_z;

   val_rot_z = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderPx(hObject,eventdata)
   global val_trans_x;

   val_trans_x = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderPy(hObject,eventdata)
   global val_trans_y;

   val_trans_y = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderPz(hObject,eventdata)
   global val_trans_z;

   val_trans_z = eventdata.Value;
   
   updateFigure();
end

function updateUiRadiobutton(hObject,eventdata)
   global val_rot_x;
   global val_rot_y;
   global val_rot_z;
   global val_trans_x;
   global val_trans_y;
   global val_trans_z;
   global scenario;
   global sld_rx;
   global sld_ry;
   global sld_rz;
   global sld_px;
   global sld_py;
   global sld_pz;
   global cbx_rx;
   global cbx_ry;
   global cbx_rz;
   global cbx_px;
   global cbx_py;
   global cbx_pz;   
   global rb_r0;
   global rb_r1;
   global rb_r2;
   global rb_p0;
   global rb_p1;
   global rb_p2;
   
   val_rot_x = 0; sld_rx.Value = val_rot_x; 
   val_rot_y = 0; sld_ry.Value = val_rot_y; 
   val_rot_z = 0; sld_rz.Value = val_rot_z; 
   val_trans_x = 0; sld_px.Value = val_trans_x; 
   val_trans_y = 0; sld_py.Value = val_trans_y; 
   val_trans_z = 0; sld_pz.Value = val_trans_z; 
           
   switch hObject.SelectedObject 
     case rb_r0
       scenario =  0;
       sld_px.Enable = 0;
       sld_py.Enable = 0;
       sld_pz.Enable = 0;
       cbx_px.Enable = 0;
       cbx_py.Enable = 0;
       cbx_pz.Enable = 0;
     case rb_r1
       scenario =  1;
       sld_px.Enable = 0;
       sld_py.Enable = 0;
       sld_pz.Enable = 0;
       cbx_px.Enable = 0;
       cbx_py.Enable = 0;
       cbx_pz.Enable = 0;
     case rb_r2
       scenario =  2;
       sld_px.Enable = 0;
       sld_py.Enable = 0;
       sld_pz.Enable = 0;
       cbx_px.Enable = 0;
       cbx_py.Enable = 0;
       cbx_pz.Enable = 0;
     case rb_p0
       scenario =  3;
       sld_px.Enable = 1;
       sld_py.Enable = 1;
       sld_pz.Enable = 1;
       cbx_px.Enable = 1;
       cbx_py.Enable = 1;
       cbx_pz.Enable = 1;
     case rb_p1
       scenario =  4;
       sld_px.Enable = 1;
       sld_py.Enable = 1;
       sld_pz.Enable = 1;
       cbx_px.Enable = 1;
       cbx_py.Enable = 1;
       cbx_pz.Enable = 1;
     case rb_p2
       scenario =  5;
       sld_px.Enable = 1;
       sld_py.Enable = 1;
       sld_pz.Enable = 1;
       cbx_px.Enable = 1;
       cbx_py.Enable = 1;
       cbx_pz.Enable = 1;
   end
   
   updateFigure();
end