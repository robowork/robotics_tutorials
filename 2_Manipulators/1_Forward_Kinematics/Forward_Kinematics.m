clear all; clc;
global val_theta_1;
global val_theta_2;
global val_theta_3;
global val_theta_4;
global val_theta_5;
global val_theta_6;
global sld_theta_1;
global sld_theta_2;
global sld_theta_3;
global sld_theta_4;
global sld_theta_5;
global sld_theta_6;
global cbx_theta_1;
global cbx_theta_2;
global cbx_theta_3;
global cbx_theta_4;
global cbx_theta_5;
global cbx_theta_6;
global rb_fk0;
global rb_fk1;
global scenario;
global fig_ax;
global h_fig_caz;
global h_fig_cel;
global S_p;
global S_R;
global S_T;
bdclose all; close all;

val_theta_1 = 0;
val_theta_2 = 0;
val_theta_3 = 0;
val_theta_4 = 0;
val_theta_5 = 0;
val_theta_6 = 0;

scenario = 0;

h_uifig = uifigure('Position',[10 10 550 850]); %,'ValueChangedFcn',@updateFigure);
sld_theta_1 = uislider('Parent',h_uifig,'Position',[100 825-0*50+0*5 400 10],'Limits',[-180 180],'Value',val_theta_1);
sld_theta_2 = uislider('Parent',h_uifig,'Position',[100 825-1*50+1*5 400 10],'Limits',[-180 180],'Value',val_theta_2);
sld_theta_3 = uislider('Parent',h_uifig,'Position',[100 825-2*50+2*5 400 10],'Limits',[-180 180],'Value',val_theta_3);
sld_theta_4 = uislider('Parent',h_uifig,'Position',[100 825-3*50+3*5 400 10],'Limits',[-180 180],'Value',val_theta_4, 'Enable',0);
sld_theta_5 = uislider('Parent',h_uifig,'Position',[100 825-4*50+4*5 400 10],'Limits',[-180 180],'Value',val_theta_5, 'Enable',0);
sld_theta_6 = uislider('Parent',h_uifig,'Position',[100 825-5*50+5*5 400 10],'Limits',[-180 180],'Value',val_theta_6, 'Enable',0);
fig_ax = uiaxes('Parent',h_uifig,'Position',[10 10 500 500]);
sld_theta_1.ValueChangingFcn = @updateUiSliderTheta1;
sld_theta_2.ValueChangingFcn = @updateUiSliderTheta2;
sld_theta_3.ValueChangingFcn = @updateUiSliderTheta3;
sld_theta_4.ValueChangingFcn = @updateUiSliderTheta4;
sld_theta_5.ValueChangingFcn = @updateUiSliderTheta5;
sld_theta_6.ValueChangingFcn = @updateUiSliderTheta6;
cbx_theta_1 = uicheckbox('Parent',h_uifig,'Position',[5 825-0*50+0*5 400 10],'Text','reset       θ1','Value',0,...
                         'ValueChangedFcn',@(cbx_theta_1,event) updateUiCheckboxTheta1(cbx_theta_1,sld_theta_1));
cbx_theta_2 = uicheckbox('Parent',h_uifig,'Position',[5 825-1*50+1*5 400 10],'Text','reset       θ2','Value',0,...
                         'ValueChangedFcn',@(cbx_theta_2,event) updateUiCheckboxTheta2(cbx_theta_2,sld_theta_2));
cbx_theta_3 = uicheckbox('Parent',h_uifig,'Position',[5 825-2*50+2*5 400 10],'Text','reset       θ3','Value',0,...
                         'ValueChangedFcn',@(cbx_theta_3,event) updateUiCheckboxTheta3(cbx_theta_3,sld_theta_3));
cbx_theta_4 = uicheckbox('Parent',h_uifig,'Position',[5 825-3*50+3*5 400 10],'Text','reset       θ4','Value',0, 'Enable',0,...
                         'ValueChangedFcn',@(cbx_theta_4,event) updateUiCheckboxTheta4(cbx_theta_4,sld_theta_4));
cbx_theta_5 = uicheckbox('Parent',h_uifig,'Position',[5 825-4*50+4*5 400 10],'Text','reset       θ5','Value',0, 'Enable',0,...
                         'ValueChangedFcn',@(cbx_theta_5,event) updateUiCheckboxTheta5(cbx_theta_5,sld_theta_5));
cbx_theta_6 = uicheckbox('Parent',h_uifig,'Position',[5 825-5*50+5*5 400 10],'Text','reset       θ6','Value',0, 'Enable',0,...
                         'ValueChangedFcn',@(cbx_theta_6,event) updateUiCheckboxTheta6(cbx_theta_6,sld_theta_6));
bg = uibuttongroup(h_uifig,'Position',[5 825-6*50+0*5 550-10 40]);
bg.SelectionChangedFcn = @updateUiRadiobutton;
rb_fk0 = uiradiobutton('Parent',bg,'Position',[5+0*75 0 100 40],'Text','FK #1');
rb_fk1 = uiradiobutton('Parent',bg,'Position',[5+1*75 0 100 40],'Text','FK #2');

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
set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-2.5 2.5],'ylim',[-2.5 2.5],'zlim',[-2.5 2.5]);
dummy_hObject.SelectedObject = rb_fk0;
updateUiRadiobutton(dummy_hObject,NaN);

[h_fig_caz,h_fig_cel] = view(fig_ax);

function updateFigure()
   global val_theta_1;
   global val_theta_2; 
   global val_theta_3;
   global val_theta_4;
   global val_theta_5;
   global val_theta_6;
   global scenario;
   global fig_ax;
   global h_fig_caz;
   global h_fig_cel;
   global S_p;
   global S_R;
   global S_T;

   [h_fig_caz,h_fig_cel] = view(fig_ax);
   
   switch scenario 
    case 0
      % 3-Joint SE(2)
      L1 = 0.75;
      L2 = 0.5;
      L3 = 0.25;

      Home_T = [eye(3)     [(L1+L2+L3);0;0];
                zeros(1,3) 1];
      
      S1 = [[0;0;1];
            cross(-[0;0;1] , [-(L1+L2+L3);0;0])]; 
      S2 = [[0;0;0];
            [1;0;0]];
      S3 = [[0;0;1];
            cross(-[0;0;1] , [-(L3);0;0])]; 

      if abs(deg2rad(val_theta_1)) < 1e-3
        e_S1Matrix_theta1 = [eye(3)     zeros(3,1);
                             zeros(1,3) 1];  
      elseif norm(S1(1:3)) < 1e-6
        V = eye(3);
        e_S1Matrix_theta1 = [eye(3)     V*S1(4:6) * deg2rad(val_theta_1);
                             zeros(1,3) 1];
      else
        omegaSkew = [0      -S1(3) S1(2);
                     S1(3)  0      -S1(1);
                     -S1(2) S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_S1Matrix_theta1 = [e_omegaSkew V*S1(4:6) * deg2rad(val_theta_1);
                             zeros(1,3)  1];
      end

      if abs(deg2rad(val_theta_2)/pi) < 1e-3
        e_S2Matrix_theta2 = [eye(3)     zeros(3,1);
                             zeros(1,3) 1];  
      elseif norm(S2(1:3)) < 1e-6
        V = eye(3);
        e_S2Matrix_theta2 = [eye(3)     V*S2(4:6) * deg2rad(val_theta_2)/pi
                             zeros(1,3) 1];
      else
        omegaSkew = [0      -S2(3) S2(2);
                     S2(3)  0      -S2(1);
                     -S2(2) S2(1)  0] * deg2rad(val_theta_2)/pi;
        theta = norm(S2(1:3) * deg2rad(val_theta_2)/pi);
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_S2Matrix_theta2 = [e_omegaSkew V*S2(4:6) * deg2rad(val_theta_2)/pi;
                             zeros(1,3)  1];
      end
      
      if abs(deg2rad(val_theta_3)) < 1e-3
        e_S3Matrix_theta3 = [eye(3)     zeros(3,1);
                             zeros(1,3) 1];  
      elseif norm(S3(1:3)) < 1e-6
        V = eye(3);
        e_S3Matrix_theta3 = [eye(3)     V*S3(4:6) * deg2rad(val_theta_3)
                             zeros(1,3) 1];
      else
        omegaSkew = [0      -S3(3) S3(2);
                     S3(3)  0      -S3(1);
                     -S3(2) S3(1)  0] * deg2rad(val_theta_3);
        theta = norm(S3(1:3) * deg2rad(val_theta_3));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_S3Matrix_theta3 = [e_omegaSkew V*S3(4:6) * deg2rad(val_theta_3);
                             zeros(1,3)  1];
      end

      E_T = ((Home_T * e_S1Matrix_theta1) * e_S2Matrix_theta2) * e_S3Matrix_theta3;  % Post-multiply

      % Intermediate frames
      % Joint 1
      Home_j1 = [eye(3)     [(0);0;0];
                 zeros(1,3) 1];
             
      j1_S1 = [[0;0;1];
               cross(-[0;0;1] , [-(0);0;0])]; 
           
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_j1S1Matrix_theta1 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j1_S1(1:3)) < 1e-6
        V = eye(3);
        e_j1S1Matrix_theta1 = [eye(3)     V*j1_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j1_S1(3) j1_S1(2);
                     j1_S1(3)  0         -j1_S1(1);
                     -j1_S1(2) j1_S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(j1_S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j1S1Matrix_theta1 = [e_omegaSkew V*j1_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3)  1];
      end
      
      j1_T = Home_j1 * e_j1S1Matrix_theta1;  % Post-multiply

      % Joint 2
      Home_j2 = [eye(3)     [(L1+L2);0;0];
                 zeros(1,3) 0];
             
      j2_S1 = [[0;0;1];
               cross(-[0;0;1] , [-(L1+L2);0;0])]; 
      j2_S2 = [[0;0;0];
               [1;0;0]];
           
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_j2S1Matrix_theta1 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j2_S1(1:3)) < 1e-6
        V = eye(3);
        e_j2S1Matrix_theta1 = [eye(3)     V*j2_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j2_S1(3) j2_S1(2);
                     j2_S1(3)  0         -j2_S1(1);
                     -j2_S1(2) j2_S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(j2_S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j2S1Matrix_theta1 = [e_omegaSkew V*j2_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3)  1];
      end

      if abs(deg2rad(val_theta_2)/pi) < 1e-3
        e_j2S2Matrix_theta2 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j2_S2(1:3)) < 1e-6
        V = eye(3);
        e_j2S2Matrix_theta2 = [eye(3)     V*j2_S2(4:6) * deg2rad(val_theta_2)/pi
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j2_S2(3) j2_S2(2);
                     j2_S2(3)  0         -j2_S2(1);
                     -j2_S2(2) j2_S2(1)  0] * deg2rad(val_theta_2)/pi;
        theta = norm(j2_S2(1:3) * deg2rad(val_theta_2)/pi);
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j2S2Matrix_theta2 = [e_omegaSkew V*j2_S2(4:6) * deg2rad(val_theta_2)/pi;
                               zeros(1,3)  1];
      end     
      
      j2_T = (Home_j2 * e_j2S1Matrix_theta1) * e_j2S2Matrix_theta2;  % Post-multiply

      % Joint 3
      Home_j3 = [eye(3)     [(L1+L2);0;0];
                 zeros(1,3) 1];
             
      j3_S1 = [[0;0;1];
               cross(-[0;0;1] , [-(L1+L2);0;0])]; 
      j3_S2 = [[0;0;0];
               [1;0;0]];
      j3_S3 = [[0;0;1];
               cross(-[0;0;1] , [-(0);0;0])]; 
        
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_j3S1Matrix_theta1 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j3_S1(1:3)) < 1e-6
        V = eye(3);
        e_j3S1Matrix_theta1 = [eye(3)     V*j3_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j3_S1(3) j3_S1(2);
                     j3_S1(3)  0         -j3_S1(1);
                     -j3_S1(2) j3_S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(j3_S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j3S1Matrix_theta1 = [e_omegaSkew V*j3_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3)  1];
      end

      if abs(deg2rad(val_theta_2)/pi) < 1e-3
        e_j3S2Matrix_theta2 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j3_S2(1:3)) < 1e-6
        V = eye(3);
        e_j3S2Matrix_theta2 = [eye(3)     V*j3_S2(4:6) * deg2rad(val_theta_2)/pi
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j3_S2(3) j3_S2(2);
                     j3_S2(3)  0         -j3_S2(1);
                     -j3_S2(2) j3_S2(1)  0] * deg2rad(val_theta_2)/pi;
        theta = norm(j3_S2(1:3) * deg2rad(val_theta_2)/pi);
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j3S2Matrix_theta2 = [e_omegaSkew V*j3_S2(4:6) * deg2rad(val_theta_2)/pi;
                               zeros(1,3)  1];
      end    
      
      if abs(deg2rad(val_theta_3)) < 1e-3
        e_j3S3Matrix_theta3 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j3_S3(1:3)) < 1e-6
        V = eye(3);
        e_j3S3Matrix_theta3 = [eye(3)     V*j3_S3(4:6) * deg2rad(val_theta_3);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j3_S3(3) j3_S3(2);
                     j3_S3(3)  0         -j3_S3(1);
                     -j3_S3(2) j3_S3(1)  0] * deg2rad(val_theta_3);
        theta = norm(j3_S3(1:3) * deg2rad(val_theta_3));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j3S3Matrix_theta3 = [e_omegaSkew V*j3_S3(4:6) * deg2rad(val_theta_3);
                               zeros(1,3)  1];
      end
      
      j3_T = ((Home_j3 * e_j3S1Matrix_theta1) * e_j3S2Matrix_theta2) * e_j3S3Matrix_theta3;  % Post-multiply

      % Result
      translations = [S_p';
                      j1_T(1:3,4)';
                      j2_T(1:3,4)';
                      j3_T(1:3,4)';
                      E_T(1:3,4)'];
      rotations = [rotm2quat(S_R);
                   rotm2quat(j1_T(1:3,1:3));
                   rotm2quat(j2_T(1:3,1:3));
                   rotm2quat(j3_T(1:3,1:3));
                   rotm2quat(E_T(1:3,1:3))];
               
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.5,'Parent',fig_ax);
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',0.75,'Parent',fig_ax);            
      plotTransforms(translations(3,:),rotations(3,:),'FrameSize',0.75,'Parent',fig_ax);            
      plotTransforms(translations(4,:),rotations(4,:),'FrameSize',0.75,'Parent',fig_ax);            
      plotTransforms(translations(5,:),rotations(5,:),'FrameSize',1.0,'Parent',fig_ax);        
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-2.5 2.5],'ylim',[-2.5 2.5],'zlim',[-2.5 2.5]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = '3-Joint SE(2) Forward Kinematics (Body Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    case 1
      % 3-Joint SE(3)
      L1 = 0.25;
      L2 = 0.10;
      L3 = 0.75;
      L4 = 0.10;
      L5 = 0.50;

      Home_T = [1 0 0   (-L4);
                0 1 0   (L2);
                0 0 1   (L1+L3+L5);
                0 0 0   1];
      
      S1 = [[0;0;1];
            cross(-[0;0;1] , [L4;-(L2);-(L1+L3+L5)])]; 
      S2 = [[0;1;0];
            cross(-[0;1;0] , [L4;-(0);-(L3+L5)])]; 
      S3 = [[1;0;0];
            cross(-[1;0;0] , [0;-(0);-(L5)])]; 
        
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_S1Matrix_theta1 = [eye(3)     zeros(3,1);
                             zeros(1,3) 1];  
      elseif norm(S1(1:3)) < 1e-6
        V = eye(3);
        e_S1Matrix_theta1 = [eye(3)     V*S1(4:6) * deg2rad(val_theta_1);
                             zeros(1,3) 1];
      else
        omegaSkew = [0      -S1(3) S1(2);
                     S1(3)  0      -S1(1);
                     -S1(2) S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_S1Matrix_theta1 = [e_omegaSkew V*S1(4:6) * deg2rad(val_theta_1);
                             zeros(1,3)  1];
      end

      if abs(deg2rad(val_theta_2)) < 1e-3
        e_S2Matrix_theta2 = [eye(3)     zeros(3,1);
                             zeros(1,3) 1];  
      elseif norm(S2(1:3)) < 1e-6
        V = eye(3);
        e_S2Matrix_theta2 = [eye(3)     V*S2(4:6) * deg2rad(val_theta_2)
                             zeros(1,3) 1];
      else
        omegaSkew = [0      -S2(3) S2(2);
                     S2(3)  0      -S2(1);
                     -S2(2) S2(1)  0] * deg2rad(val_theta_2);
        theta = norm(S2(1:3) * deg2rad(val_theta_2));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_S2Matrix_theta2 = [e_omegaSkew V*S2(4:6) * deg2rad(val_theta_2);
                             zeros(1,3)  1];
      end
      
      if abs(deg2rad(val_theta_3)) < 1e-3
        e_S3Matrix_theta3 = [eye(3)     zeros(3,1);
                             zeros(1,3) 1];  
      elseif norm(S3(1:3)) < 1e-6
        V = eye(3);
        e_S3Matrix_theta3 = [eye(3)     V*S3(4:6) * deg2rad(val_theta_3)
                             zeros(1,3) 1];
      else
        omegaSkew = [0      -S3(3) S3(2);
                     S3(3)  0      -S3(1);
                     -S3(2) S3(1)  0] * deg2rad(val_theta_3);
        theta = norm(S3(1:3) * deg2rad(val_theta_3));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_S3Matrix_theta3 = [e_omegaSkew V*S3(4:6) * deg2rad(val_theta_3);
                             zeros(1,3)  1];
      end

      E_T = ((Home_T * e_S1Matrix_theta1) * e_S2Matrix_theta2) * e_S3Matrix_theta3;  % Post-multiply
      
      % Intermediate frames
      % Joint 1
      Home_j1 = [1 0 0   (-0);
                 0 1 0   (0);
                 0 0 1   (0);
                 0 0 0   1];
             
      j1_S1 = [[0;0;1];
               cross(-[0;0;1] , [-(0);0;0])]; 
           
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_j1S1Matrix_theta1 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j1_S1(1:3)) < 1e-6
        V = eye(3);
        e_j1S1Matrix_theta1 = [eye(3)     V*j1_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j1_S1(3) j1_S1(2);
                     j1_S1(3)  0         -j1_S1(1);
                     -j1_S1(2) j1_S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(j1_S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j1S1Matrix_theta1 = [e_omegaSkew V*j1_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3)  1];
      end
      
      j1_T = Home_j1 * e_j1S1Matrix_theta1;  % Post-multiply
      
      % Joint 2
      Home_j2 = [1 0 0   (-0);
                 0 1 0   (L2);
                 0 0 1   (L1);
                 0 0 0   1];
             
      j2_S1 = [[0;0;1];
               cross(-[0;0;1] , [0;-(L2);-(L1)])]; 
      j2_S2 = [[0;1;0];
               cross(-[0;1;0] , [0;-(0);-(0)])]; 
           
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_j2S1Matrix_theta1 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j2_S1(1:3)) < 1e-6
        V = eye(3);
        e_j2S1Matrix_theta1 = [eye(3)     V*j2_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j2_S1(3) j2_S1(2);
                     j2_S1(3)  0         -j2_S1(1);
                     -j2_S1(2) j2_S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(j2_S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j2S1Matrix_theta1 = [e_omegaSkew V*j2_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3)  1];
      end

      if abs(deg2rad(val_theta_2)) < 1e-3
        e_j2S2Matrix_theta2 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j2_S2(1:3)) < 1e-6
        V = eye(3);
        e_j2S2Matrix_theta2 = [eye(3)     V*j2_S2(4:6) * deg2rad(val_theta_2)
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j2_S2(3) j2_S2(2);
                     j2_S2(3)  0         -j2_S2(1);
                     -j2_S2(2) j2_S2(1)  0] * deg2rad(val_theta_2);
        theta = norm(j2_S2(1:3) * deg2rad(val_theta_2));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j2S2Matrix_theta2 = [e_omegaSkew V*j2_S2(4:6) * deg2rad(val_theta_2);
                               zeros(1,3)  1];
      end     
      
      j2_T = (Home_j2 * e_j2S1Matrix_theta1) * e_j2S2Matrix_theta2;  % Post-multiply
      
      % Joint 3
      Home_j3 = [1 0 0   (-L4);
                 0 1 0   (L2);
                 0 0 1   (L1+L3);
                 0 0 0   1];
             
      j3_S1 = [[0;0;1];
               cross(-[0;0;1] , [L4;-(L2);-(L1+L3)])]; 
      j3_S2 = [[0;1;0];
               cross(-[0;1;0] , [L4;-(0);-(L3)])]; 
      j3_S3 = [[1;0;0];
               cross(-[1;0;0] , [0;-(0);-(0)])]; 
        
      if abs(deg2rad(val_theta_1)) < 1e-3
        e_j3S1Matrix_theta1 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j3_S1(1:3)) < 1e-6
        V = eye(3);
        e_j3S1Matrix_theta1 = [eye(3)     V*j3_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j3_S1(3) j3_S1(2);
                     j3_S1(3)  0         -j3_S1(1);
                     -j3_S1(2) j3_S1(1)  0] * deg2rad(val_theta_1);
        theta = norm(j3_S1(1:3) * deg2rad(val_theta_1));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j3S1Matrix_theta1 = [e_omegaSkew V*j3_S1(4:6) * deg2rad(val_theta_1);
                               zeros(1,3)  1];
      end

      if abs(deg2rad(val_theta_2)) < 1e-3
        e_j3S2Matrix_theta2 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j3_S2(1:3)) < 1e-6
        V = eye(3);
        e_j3S2Matrix_theta2 = [eye(3)     V*j3_S2(4:6) * deg2rad(val_theta_2)
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j3_S2(3) j3_S2(2);
                     j3_S2(3)  0         -j3_S2(1);
                     -j3_S2(2) j3_S2(1)  0] * deg2rad(val_theta_2);
        theta = norm(j3_S2(1:3) * deg2rad(val_theta_2));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j3S2Matrix_theta2 = [e_omegaSkew V*j3_S2(4:6) * deg2rad(val_theta_2);
                               zeros(1,3)  1];
      end    
      
      if abs(deg2rad(val_theta_3)) < 1e-3
        e_j3S3Matrix_theta3 = [eye(3)     zeros(3,1);
                               zeros(1,3) 1];  
      elseif norm(j3_S3(1:3)) < 1e-6
        V = eye(3);
        e_j3S3Matrix_theta3 = [eye(3)     V*j3_S3(4:6) * deg2rad(val_theta_3);
                               zeros(1,3) 1];
      else
        omegaSkew = [0         -j3_S3(3) j3_S3(2);
                     j3_S3(3)  0         -j3_S3(1);
                     -j3_S3(2) j3_S3(1)  0] * deg2rad(val_theta_3);
        theta = norm(j3_S3(1:3) * deg2rad(val_theta_3));
        e_omegaSkew = eye(3) + ( sin(theta)/theta ) * omegaSkew + ( (1-cos(theta))/(theta^2) ) * omegaSkew^2;
        V = eye(3) + ( (1-cos(theta))/(theta^2) ) * omegaSkew + ( (theta-sin(theta))/(theta^3) ) * omegaSkew^2;
        e_j3S3Matrix_theta3 = [e_omegaSkew V*j3_S3(4:6) * deg2rad(val_theta_3);
                               zeros(1,3)  1];
      end
      
      j3_T = ((Home_j3 * e_j3S1Matrix_theta1) * e_j3S2Matrix_theta2) * e_j3S3Matrix_theta3;  % Post-multiply
      
      % Result
      translations = [S_p';
                      j1_T(1:3,4)';
                      j2_T(1:3,4)';
                      j3_T(1:3,4)';
                      E_T(1:3,4)'];
      rotations = [rotm2quat(S_R);
                   rotm2quat(j1_T(1:3,1:3));
                   rotm2quat(j2_T(1:3,1:3));
                   rotm2quat(j3_T(1:3,1:3));
                   rotm2quat(E_T(1:3,1:3));];
               
      cla(fig_ax);
      plotTransforms(translations(1,:),rotations(1,:),'FrameSize',0.5,'Parent',fig_ax);          
      plotTransforms(translations(2,:),rotations(2,:),'FrameSize',0.75,'Parent',fig_ax);          
      plotTransforms(translations(3,:),rotations(3,:),'FrameSize',0.75,'Parent',fig_ax);          
      plotTransforms(translations(4,:),rotations(4,:),'FrameSize',0.75,'Parent',fig_ax);          
      plotTransforms(translations(5,:),rotations(5,:),'FrameSize',1.0,'Parent',fig_ax);      
      set(fig_ax,'dataaspectratio',[1 1 1],'xgrid',1,'ygrid',1,'zgrid',1,'xlim',[-2.5 2.5],'ylim',[-2.5 2.5],'zlim',[-2.5 2.5]);
      fig_ax_text = get(fig_ax,'title'); 
      fig_ax_text.String = '2-Joint SE(3) Forward Kinematics (Body Frame)';
      set(fig_ax,'title',fig_ax_text); 
      view(fig_ax,h_fig_caz,h_fig_cel);
      
    otherwise  
      bdclose all; close all;
 
   end

end

function updateUiCheckboxTheta1(hObject,sldObject)
   global val_theta_1;
   
   if hObject.Value
       val_theta_1 = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxTheta2(hObject,sldObject)
   global val_theta_2;
   
   if hObject.Value
       val_theta_2 = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxTheta3(hObject,sldObject)
   global val_theta_3;
   
   if hObject.Value
       val_theta_3 = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxTheta4(hObject,sldObject)
   global val_theta_4;
   
   if hObject.Value
       val_theta_4 = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxTheta5(hObject,sldObject)
   global val_theta_5;
   
   if hObject.Value
       val_theta_5 = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end
function updateUiCheckboxTheta6(hObject,sldObject)
   global val_theta_6;
   
   if hObject.Value
       val_theta_6 = 0;

       updateFigure();

       hObject.Value = 0;
       sldObject.Value = 0;
   end
end

function updateUiSliderTheta1(hObject,eventdata)
   global val_theta_1;

   val_theta_1 = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderTheta2(hObject,eventdata)
   global val_theta_2;

   val_theta_2 = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderTheta3(hObject,eventdata)
   global val_theta_3;

   val_theta_3 = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderTheta4(hObject,eventdata)
   global val_theta_4;

   val_theta_4 = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderTheta5(hObject,eventdata)
   global val_theta_5;

   val_theta_5 = eventdata.Value;
   
   updateFigure();
end
function updateUiSliderTheta6(hObject,eventdata)
   global val_theta_6;

   val_theta_6 = eventdata.Value;
   
   updateFigure();
end

function updateUiRadiobutton(hObject,eventdata)
   global val_theta_1;
   global val_theta_2;
   global val_theta_3;
   global val_theta_4;
   global val_theta_5;
   global val_theta_6;
   global sld_theta_1;
   global sld_theta_2;
   global sld_theta_3;
   global sld_theta_4;
   global sld_theta_5;
   global sld_theta_6;
   global cbx_theta_1;
   global cbx_theta_2;
   global cbx_theta_3;
   global cbx_theta_4;
   global cbx_theta_5;
   global cbx_theta_6;
   global scenario;
   global rb_fk0;
   global rb_fk1;
   
   val_theta_1 = 0; sld_theta_1.Value = val_theta_1; 
   val_theta_2 = 0; sld_theta_2.Value = val_theta_2; 
   val_theta_3 = 0; sld_theta_3.Value = val_theta_3; 
   val_theta_4 = 0; sld_theta_4.Value = val_theta_4; 
   val_theta_5 = 0; sld_theta_5.Value = val_theta_5; 
   val_theta_6 = 0; sld_theta_6.Value = val_theta_6; 
           
   switch hObject.SelectedObject 
     case rb_fk0
       scenario =  0;
       sld_theta_1.Enable = 1;
       sld_theta_2.Enable = 1;
       sld_theta_3.Enable = 1;
       sld_theta_4.Enable = 0;
       sld_theta_5.Enable = 0;
       sld_theta_6.Enable = 0;
       cbx_theta_1.Enable = 1;
       cbx_theta_2.Enable = 1;
       cbx_theta_3.Enable = 1;
       cbx_theta_4.Enable = 0;
       cbx_theta_5.Enable = 0;
       cbx_theta_6.Enable = 0;
     case rb_fk1
       scenario =  1;
       sld_theta_1.Enable = 1;
       sld_theta_2.Enable = 1;
       sld_theta_3.Enable = 1;
       sld_theta_4.Enable = 0;
       sld_theta_5.Enable = 0;
       sld_theta_6.Enable = 0;
       cbx_theta_1.Enable = 1;
       cbx_theta_2.Enable = 1;
       cbx_theta_3.Enable = 1;
       cbx_theta_4.Enable = 0;
       cbx_theta_5.Enable = 0;
       cbx_theta_6.Enable = 0;
   end
   
   updateFigure();
end