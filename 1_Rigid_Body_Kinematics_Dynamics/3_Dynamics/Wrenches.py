import sys, os, random
import numpy as np

from transforms3d import *  #supersedes deprecated transformations.py

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import cm

from mpl_toolkits import mplot3d


class DoubleSlider(QSlider):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.decimals = 3
        self._max_int = 10 ** self.decimals

        super().setMinimum(0)
        super().setMaximum(self._max_int)

        self._min_value = 0.0
        self._max_value = 1.0

    @property
    def _value_range(self):
        return self._max_value - self._min_value

    def value(self):
        return float(super().value()) / self._max_int * self._value_range + self._min_value

    def setValue(self, value):
        super().setValue(int((value - self._min_value) / self._value_range * self._max_int))

    def setMinimum(self, value):
        if value > self._max_value:
            raise ValueError("Minimum limit cannot be higher than maximum")

        self._min_value = value
        self.setValue(self.value())

    def setMaximum(self, value):
        if value < self._min_value:
            raise ValueError("Minimum limit cannot be higher than maximum")

        self._max_value = value
        self.setValue(self.value())

    def minimum(self):
        return self._min_value

    def maximum(self):
        return self._max_value


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Rigid Body Motion: Wrenches & Dynamics')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()


        self.val_f_x = 0
        self.val_f_y = 0
        self.val_f_z = 0
        self.val_m_x = 0
        self.val_m_y = 0
        self.val_m_z = 0

        self.B_p_x_0 = 0.50
        self.B_p_y_0 = 0.25
        self.B_p_z_0 = 0.75
        self.B_p_qw_0 = 0.872
        self.B_p_qx_0 = 0.215
        self.B_p_qy_0 = 0.189
        self.B_p_qz_0 = 0.398

        self.ef_tx.setText(str(self.B_p_x_0))
        self.ef_ty.setText(str(self.B_p_y_0))
        self.ef_tz.setText(str(self.B_p_z_0))
        self.ef_qw.setText(str(self.B_p_qw_0))
        self.ef_qx.setText(str(self.B_p_qx_0))
        self.ef_qy.setText(str(self.B_p_qy_0))
        self.ef_qz.setText(str(self.B_p_qz_0))

        self.timer_period = 0.10
        self.timer_value = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_timer)
        self.timer.start(self.timer_period * 1000)

        self.omega_x_0 = 0*30.0  #deg/s
        self.omega_y_0 = 0*15.0  #deg/s
        self.omega_z_0 = 0*90.0  #deg/s
        self.vel_x_0 = 0*-0.5
        self.vel_y_0 = 0*-2.0
        self.vel_z_0 = 0*-0.25

        self.B_p_x = self.B_p_x_0
        self.B_p_y = self.B_p_y_0
        self.B_p_z = self.B_p_z_0
        self.B_p_qw = self.B_p_qw_0  
        self.B_p_qx = self.B_p_qx_0
        self.B_p_qy = self.B_p_qy_0
        self.B_p_qz = self.B_p_qz_0 
        self.omega_x = self.omega_x_0  
        self.omega_y = self.omega_y_0
        self.omega_z = self.omega_z_0
        self.vel_x = self.vel_x_0
        self.vel_y = self.vel_y_0
        self.vel_z = self.vel_z_0
        self.dot_omega_x = 0  
        self.dot_omega_y = 0
        self.dot_omega_z = 0
        self.dot_vel_x = 0
        self.dot_vel_y = 0
        self.dot_vel_z = 0

        self.scenario = 0
        
        # rigid body component (emulated as ellipsoids) masses and semiaxes
        ellipsoid_1_mass = 1.0
        ellipsoid_1_a = 0.25
        ellipsoid_1_b = 0.75
        ellipsoid_1_c = 0.15
        ellipsoid_1_angle_z = np.deg2rad(45)
        ellipsoid_1_angle_y = np.deg2rad(30)
        ellipsoid_1_angle_x = np.deg2rad(15)

        # Spherical angles
        ellipsoid_1_u = np.linspace(0, 2*np.pi, 25)
        ellipsoid_1_v = np.linspace(0, np.pi, 25)
        # Cartesian coordinates
        self.ellipsoid_1_x = ellipsoid_1_a * np.outer(np.cos(ellipsoid_1_u), np.sin(ellipsoid_1_v))
        self.ellipsoid_1_y = ellipsoid_1_b * np.outer(np.sin(ellipsoid_1_u), np.sin(ellipsoid_1_v))
        self.ellipsoid_1_z = ellipsoid_1_c * np.outer(np.ones_like(ellipsoid_1_u), np.cos(ellipsoid_1_v))

        ellipsoid_1_I = (ellipsoid_1_mass/5)*np.array([[ellipsoid_1_b**2+ellipsoid_1_c**2, 0, 0], [0, ellipsoid_1_a**2+ellipsoid_1_c**2, 0], [0, 0, ellipsoid_1_a**2+ellipsoid_1_b**2]])
        ellipsoid_1_R_z = np.array([[np.cos(ellipsoid_1_angle_z), -np.sin(ellipsoid_1_angle_z), 0], [np.sin(ellipsoid_1_angle_z), np.cos(ellipsoid_1_angle_z), 0], [0, 0, 1]])
        ellipsoid_1_R_y = np.array([[np.cos(ellipsoid_1_angle_y), 0, np.sin(ellipsoid_1_angle_y)], [0, 1, 0], [-np.sin(ellipsoid_1_angle_y), 0, np.cos(ellipsoid_1_angle_y)]])
        ellipsoid_1_R_x = np.array([[1, 0, 0], [0, np.cos(ellipsoid_1_angle_x), -np.sin(ellipsoid_1_angle_x)], [0, np.sin(ellipsoid_1_angle_x), np.cos(ellipsoid_1_angle_x)]])
        ellipsoid_1_R_zyx = (ellipsoid_1_R_z.dot(ellipsoid_1_R_y)).dot(ellipsoid_1_R_x)
        ellipsoid_1_I = ellipsoid_1_R_zyx * ellipsoid_1_I * np.transpose(ellipsoid_1_R_zyx)
        
        ellipsoid_2_mass = 5.0
        ellipsoid_2_a = 1.0
        ellipsoid_2_b = 0.10
        ellipsoid_2_c = 0.50
        ellipsoid_2_angle_z = np.deg2rad(-30)
        ellipsoid_2_angle_y = np.deg2rad(-15)
        ellipsoid_2_angle_x = np.deg2rad(45)
        
        # Spherical angles
        ellipsoid_2_u = np.linspace(0, 2*np.pi, 25)
        ellipsoid_2_v = np.linspace(0, np.pi, 25)
        # Cartesian coordinates
        self.ellipsoid_2_x = ellipsoid_2_a * np.outer(np.cos(ellipsoid_2_u), np.sin(ellipsoid_2_v))
        self.ellipsoid_2_y = ellipsoid_2_b * np.outer(np.sin(ellipsoid_2_u), np.sin(ellipsoid_2_v))
        self.ellipsoid_2_z = ellipsoid_2_c * np.outer(np.ones_like(ellipsoid_2_u), np.cos(ellipsoid_2_v))

        ellipsoid_2_I = (ellipsoid_2_mass/5)*np.array([[ellipsoid_2_b**2+ellipsoid_2_c**2, 0, 0], [0, ellipsoid_2_a**2+ellipsoid_2_c**2, 0], [0, 0, ellipsoid_2_a**2+ellipsoid_2_b**2]])
        ellipsoid_2_R_z = np.array([[np.cos(ellipsoid_2_angle_z), -np.sin(ellipsoid_2_angle_z), 0], [np.sin(ellipsoid_2_angle_z), np.cos(ellipsoid_2_angle_z), 0], [0, 0, 1]])
        ellipsoid_2_R_y = np.array([[np.cos(ellipsoid_2_angle_y), 0, np.sin(ellipsoid_2_angle_y)], [0, 1, 0], [-np.sin(ellipsoid_2_angle_y), 0, np.cos(ellipsoid_2_angle_y)]])
        ellipsoid_2_R_x = np.array([[1, 0, 0], [0, np.cos(ellipsoid_2_angle_x), -np.sin(ellipsoid_2_angle_x)], [0, np.sin(ellipsoid_2_angle_x), np.cos(ellipsoid_2_angle_x)]])
        ellipsoid_2_R_zyx = (ellipsoid_2_R_z.dot(ellipsoid_2_R_y)).dot(ellipsoid_2_R_x)
        ellipsoid_2_I = ellipsoid_2_R_zyx * ellipsoid_2_I * np.transpose(ellipsoid_2_R_zyx)

        # combined rigid body mass and inertia tensor
        self.combined_mass = ellipsoid_1_mass + ellipsoid_2_mass
        self.combined_I = ellipsoid_1_I + ellipsoid_2_I

        # Space Frame / origin
        self.S_p = np.array([[0], \
                             [0], \
                             [0]])
        self.S_R = np.array([[1, 0, 0], \
                             [0, 1, 0], \
                             [0, 0, 1]])
        self.S_T = np.concatenate((np.concatenate((self.S_R, self.S_p), axis=1), \
                                   np.array([[0, 0, 0, 1]])), axis=0)
   
        self.quiver_Sx = np.array([[1], \
                                   [0], \
                                   [0]])
        self.quiver_Sy = np.array([[0], \
                                   [1], \
                                   [0]])
        self.quiver_Sz = np.array([[0], \
                                   [0], \
                                   [1]])
        self.translation_S = self.S_p
        self.rotation_S = quaternions.mat2quat(self.S_R)

        #self.on_draw()

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_timer(self): 
        self.timer_value = self.timer_value + self.timer_period
        if self.timer_value > 5:
            self.timer_value = 0

            self.B_p_x = self.B_p_x_0
            self.B_p_y = self.B_p_y_0
            self.B_p_z = self.B_p_z_0
            self.B_p_qw = self.B_p_qw_0  
            self.B_p_qx = self.B_p_qx_0
            self.B_p_qy = self.B_p_qy_0
            self.B_p_qz = self.B_p_qz_0 
            self.omega_x = self.omega_x_0  
            self.omega_y = self.omega_y_0
            self.omega_z = self.omega_z_0
            self.vel_x = self.vel_x_0
            self.vel_y = self.vel_y_0
            self.vel_z = self.vel_z_0
            self.dot_omega_x = 0  
            self.dot_omega_y = 0
            self.dot_omega_z = 0
            self.dot_vel_x = 0
            self.dot_vel_y = 0
            self.dot_vel_z = 0
        #print(self.timer_value)

        self.on_draw()

    def on_draw(self):
        self.axes.clear()        
        self.axes.grid(True)

        if self.scenario == 0:
            # Twist expressed w.r.t. Body Frame
            B_p_zyx = np.array([[self.B_p_x], \
                                [self.B_p_y], \
                                [self.B_p_z]])
            B_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), B_p_zyx), axis=1), \
                                            np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx_rot = np.concatenate((np.concatenate((quaternions.quat2mat(np.array([self.B_p_qw, self.B_p_qx, self.B_p_qy, self.B_p_qz])), np.zeros((3,1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx = (self.S_T.dot(B_T_zyx_trans)).dot(B_T_zyx_rot)

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            ### Current Twist - Use to update KINEMATIC state ###
            Twist_B = np.array([[np.deg2rad(self.omega_x)], \
                                [np.deg2rad(self.omega_y)], \
                                [np.deg2rad(self.omega_z)], \
                                [self.vel_x], \
                                [self.vel_y], \
                                [self.vel_z]])

            if np.linalg.norm(Twist_B[0:3]) < 1e-6:
                V = np.eye(3)
                e_VMatrix_theta = np.concatenate((np.concatenate((np.eye(3), V.dot(Twist_B[3:6]) * self.timer_period), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[             0, -Twist_B[2][0],  Twist_B[1][0]], \
                                      [ Twist_B[2][0],              0, -Twist_B[0][0]], \
                                      [-Twist_B[1][0],  Twist_B[0][0],              0]]) * self.timer_period
                theta = np.linalg.norm(Twist_B[0:3] * self.timer_period)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_VMatrix_theta = np.concatenate((np.concatenate((e_omegaSkew, V.dot(Twist_B[3:6] * self.timer_period)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
      
            # Forward kinematics - compute updated pose
            B_T = B_T_zyx.dot(e_VMatrix_theta)  # Post-multiply

            self.B_p_x = B_T[0][3]
            self.B_p_y = B_T[1][3]
            self.B_p_z = B_T[2][3]
            B_p_q = quaternions.mat2quat(B_T[0:3,0:3])
            self.B_p_qw = B_p_q[0]
            self.B_p_qx = B_p_q[1]
            self.B_p_qy = B_p_q[2]
            self.B_p_qz = B_p_q[3]

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            ### Current Wrench - Use to update DYNAMIC state ###
            Wrench_B = np.array([[self.val_m_x], \
                                 [self.val_m_y], \
                                 [self.val_m_z], \
                                 [self.val_f_x], \
                                 [self.val_f_y], \
                                 [self.val_f_z]])

            # 6x6 Inertia
            G = np.concatenate((np.concatenate((self.combined_I, np.zeros((3,3))), axis=1), \
                                np.concatenate((np.zeros((3,3)), self.combined_mass * np.eye(3)), axis=1)), axis=0)
            # Twist
            omega_B_Skew = np.array([[             0, -Twist_B[2][0],  Twist_B[1][0]], \
                                     [ Twist_B[2][0],              0, -Twist_B[0][0]], \
                                     [-Twist_B[1][0],  Twist_B[0][0],              0]])
            v_B_Skew = np.array([[             0, -Twist_B[5][0],  Twist_B[4][0]], \
                                 [ Twist_B[5][0],              0, -Twist_B[3][0]], \
                                 [-Twist_B[4][0],  Twist_B[3][0],              0]])
            adv_B_transp = np.concatenate((np.concatenate((omega_B_Skew, v_B_Skew), axis=1), \
                                           np.concatenate((np.zeros((3,3)), omega_B_Skew), axis=1)), axis=0)
      
            # Forward Dynamics - compute updated twist
            dot_Twist_B = np.linalg.solve(G, Wrench_B + adv_B_transp.dot(G).dot(Twist_B))
      
            Twist_B = Twist_B + dot_Twist_B * self.timer_period
            self.omega_x = np.rad2deg(Twist_B[0][0])
            self.omega_y = np.rad2deg(Twist_B[1][0])
            self.omega_z = np.rad2deg(Twist_B[2][0])
            self.vel_x = Twist_B[3][0]
            self.vel_y = Twist_B[4][0]
            self.vel_z = Twist_B[5][0]

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bx_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bx_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bx_zyx[2]-self.quiver_Bp_zyx[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_By_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_By_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_By_zyx[2]-self.quiver_Bp_zyx[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bz_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bz_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bz_zyx[2]-self.quiver_Bp_zyx[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-3.0, 3.0)
            self.axes.set_ylim(-3.0, 3.0)
            self.axes.set_zlim(-3.0, 3.0)

            ellipsoid_1_x_transformed = self.ellipsoid_1_x.copy()
            ellipsoid_1_y_transformed = self.ellipsoid_1_y.copy()
            ellipsoid_1_z_transformed = self.ellipsoid_1_z.copy()
            for j in range(len(ellipsoid_1_x_transformed)):
                for i in range(len(ellipsoid_1_x_transformed[j])):
                    ellipsoid_1_i_xyz = B_T.dot(np.array([[ellipsoid_1_x_transformed[i][j]], [ellipsoid_1_y_transformed[i][j]], [ellipsoid_1_z_transformed[i][j]], [1]]))
                    ellipsoid_1_x_transformed[i][j] = ellipsoid_1_i_xyz[0][0]
                    ellipsoid_1_y_transformed[i][j] = ellipsoid_1_i_xyz[1][0]
                    ellipsoid_1_z_transformed[i][j] = ellipsoid_1_i_xyz[2][0]
            self.axes.plot_surface(ellipsoid_1_x_transformed, ellipsoid_1_y_transformed, ellipsoid_1_z_transformed, rstride=4, cstride=4, cmap=cm.jet)

            ellipsoid_2_x_transformed = self.ellipsoid_2_x.copy()
            ellipsoid_2_y_transformed = self.ellipsoid_2_y.copy()
            ellipsoid_2_z_transformed = self.ellipsoid_2_z.copy()
            for j in range(len(ellipsoid_2_x_transformed)):
                for i in range(len(ellipsoid_2_x_transformed[j])):
                    ellipsoid_2_i_xyz = B_T.dot(np.array([[ellipsoid_2_x_transformed[i][j]], [ellipsoid_2_y_transformed[i][j]], [ellipsoid_2_z_transformed[i][j]], [1]]))
                    ellipsoid_2_x_transformed[i][j] = ellipsoid_2_i_xyz[0][0]
                    ellipsoid_2_y_transformed[i][j] = ellipsoid_2_i_xyz[1][0]
                    ellipsoid_2_z_transformed[i][j] = ellipsoid_2_i_xyz[2][0]
            self.axes.plot_surface(ellipsoid_2_x_transformed, ellipsoid_2_y_transformed, ellipsoid_2_z_transformed, rstride=4, cstride=4, cmap=cm.jet)

        self.canvas.draw()

    def on_update_values(self):       
        self.B_p_x_0 = float(self.ef_tx.text())
        self.B_p_y_0 = float(self.ef_ty.text())
        self.B_p_z_0 = float(self.ef_tz.text())
        self.B_p_qw_0 = float(self.ef_qw.text())
        self.B_p_qx_0 = float(self.ef_qx.text())
        self.B_p_qy_0 = float(self.ef_qy.text())
        self.B_p_qz_0 = float(self.ef_qz.text())

        if self.cbx_fx.isChecked():
            self.cbx_fx.setChecked(False)
            self.sld_fx.setValue(0.0)
        if self.cbx_fy.isChecked():
            self.cbx_fy.setChecked(False)
            self.sld_fy.setValue(0.0)
        if self.cbx_fz.isChecked():
            self.cbx_fz.setChecked(False)
            self.sld_fz.setValue(0.0) 
        if self.cbx_mx.isChecked():
            self.cbx_mx.setChecked(False)
            self.sld_mx.setValue(0.0)
        if self.cbx_my.isChecked():
            self.cbx_my.setChecked(False)
            self.sld_my.setValue(0.0)
        if self.cbx_mz.isChecked():
            self.cbx_mz.setChecked(False)
            self.sld_mz.setValue(0.0)

        self.val_f_x = self.sld_fx.value()
        self.val_f_y = self.sld_fy.value()
        self.val_f_z = self.sld_fz.value()
        self.val_m_x = self.sld_mx.value()
        self.val_m_y = self.sld_my.value()
        self.val_m_z = self.sld_mz.value()

        #if self.rb_t0.isChecked() and self.scenario != 0:
        #    self.scenario = 0

        #self.on_draw()

    def create_main_frame(self):
        self.main_frame = QWidget()

        self.dpi = 100
        self.fig = Figure((5.0, 15.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        self.axes = self.fig.add_subplot(111, projection='3d', proj_type='ortho')
        self.axes.set_aspect('equal')
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        

        self.ef_tx = QLineEdit()
        self.ef_tx.setMinimumWidth(50)
        self.ef_tx.setFixedWidth(50)
        self.connect(self.ef_tx, SIGNAL('editingFinished()'), self.on_update_values)
        
        self.ef_ty = QLineEdit()
        self.ef_ty.setMinimumWidth(50)
        self.ef_ty.setFixedWidth(50)
        self.connect(self.ef_ty, SIGNAL('editingFinished()'), self.on_update_values)

        self.ef_tz = QLineEdit()
        self.ef_tz.setMinimumWidth(50)
        self.ef_tz.setFixedWidth(50)
        self.connect(self.ef_tz, SIGNAL('editingFinished()'), self.on_update_values)

        self.ef_qw = QLineEdit()
        self.ef_qw.setMinimumWidth(50)
        self.ef_qw.setFixedWidth(50)
        self.connect(self.ef_qw, SIGNAL('editingFinished()'), self.on_update_values)

        self.ef_qx = QLineEdit()
        self.ef_qx.setMinimumWidth(50)
        self.ef_qx.setFixedWidth(50)
        self.connect(self.ef_qx, SIGNAL('editingFinished()'), self.on_update_values)

        self.ef_qy = QLineEdit()
        self.ef_qy.setMinimumWidth(50)
        self.ef_qy.setFixedWidth(50)
        self.connect(self.ef_qy, SIGNAL('editingFinished()'), self.on_update_values)

        self.ef_qz = QLineEdit()
        self.ef_qz.setMinimumWidth(50)
        self.ef_qz.setFixedWidth(50)
        self.connect(self.ef_qz, SIGNAL('editingFinished()'), self.on_update_values)

        #self.draw_button = QPushButton("&Draw")
        #self.connect(self.draw_button, SIGNAL('clicked()'), self.on_update_values)
        
        self.cbx_fx = QCheckBox('reset')
        self.cbx_fx.setChecked(False)
        self.connect(self.cbx_fx, SIGNAL('stateChanged(int)'), self.on_update_values)
        
        self.sld_fx = DoubleSlider(Qt.Horizontal)
        self.sld_fx.setMinimum(-10.0)
        self.sld_fx.setMaximum(10.0)
        self.sld_fx.setValue(0.0)
        self.sld_fx.setTracking(True)
        self.sld_fx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_fx, SIGNAL('valueChanged(int)'), self.on_update_values)
        
        self.cbx_fy = QCheckBox('reset')
        self.cbx_fy.setChecked(False)
        self.connect(self.cbx_fy, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_fy = DoubleSlider(Qt.Horizontal)
        self.sld_fy.setMinimum(-10.0)
        self.sld_fy.setMaximum(10.0)
        self.sld_fy.setValue(0.0)
        self.sld_fy.setTracking(True)
        self.sld_fy.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_fy, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_fz = QCheckBox('reset')
        self.cbx_fz.setChecked(False)
        self.connect(self.cbx_fz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_fz = DoubleSlider(Qt.Horizontal)
        self.sld_fz.setMinimum(-10.0)
        self.sld_fz.setMaximum(10.0)
        self.sld_fz.setValue(0.0)
        self.sld_fz.setTracking(True)
        self.sld_fz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_fz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_mx = QCheckBox('reset')
        self.cbx_mx.setChecked(False)
        self.connect(self.cbx_mx, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_mx = DoubleSlider(Qt.Horizontal)
        self.sld_mx.setMinimum(-1.0)
        self.sld_mx.setMaximum(1.0)
        self.sld_mx.setValue(0.0)
        self.sld_mx.setTracking(True)
        self.sld_mx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_mx, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_my = QCheckBox('reset')
        self.cbx_my.setChecked(False)
        self.connect(self.cbx_my, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_my = DoubleSlider(Qt.Horizontal)
        self.sld_my.setMinimum(-1.0)
        self.sld_my.setMaximum(1.0)
        self.sld_my.setValue(0.0)
        self.sld_my.setTracking(True)
        self.sld_my.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_my, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_mz = QCheckBox('reset')
        self.cbx_mz.setChecked(False)
        self.connect(self.cbx_mz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_mz = DoubleSlider(Qt.Horizontal)
        self.sld_mz.setMinimum(-1.0)
        self.sld_mz.setMaximum(1.0)
        self.sld_mz.setValue(0.0)
        self.sld_mz.setTracking(True)
        self.sld_mz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_mz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.rb_t0 = QCheckBox('Wrench #1')
        self.rb_t0.setChecked(True)
        self.connect(self.rb_t0, SIGNAL('stateChanged(int)'), self.on_update_values)

        hbox_Rx = QHBoxLayout()
        for w in [ self.cbx_fx, QLabel('fx'), QLabel('-10'), self.sld_fx, QLabel('10')]:
            hbox_Rx.addWidget(w)
            hbox_Rx.setAlignment(w, Qt.AlignVCenter)

        hbox_Ry = QHBoxLayout()
        for w in [ self.cbx_fy, QLabel('fy'), QLabel('-10'), self.sld_fy, QLabel('10')]:
            hbox_Ry.addWidget(w)
            hbox_Ry.setAlignment(w, Qt.AlignVCenter)

        hbox_Rz = QHBoxLayout()
        for w in [ self.cbx_fz, QLabel('fz'), QLabel('-10'), self.sld_fz, QLabel('10')]:
            hbox_Rz.addWidget(w)
            hbox_Rz.setAlignment(w, Qt.AlignVCenter)

        hbox_px = QHBoxLayout()
        for w in [ self.cbx_mx, QLabel('mx'), QLabel('-1'), self.sld_mx, QLabel('1')]:
            hbox_px.addWidget(w)
            hbox_px.setAlignment(w, Qt.AlignVCenter)

        hbox_py = QHBoxLayout()
        for w in [ self.cbx_my, QLabel('my'), QLabel('-1'), self.sld_my, QLabel('1')]:
            hbox_py.addWidget(w)
            hbox_py.setAlignment(w, Qt.AlignVCenter)

        hbox_pz = QHBoxLayout()
        for w in [ self.cbx_mz, QLabel('mz'), QLabel('-1'), self.sld_mz, QLabel('1')]:
            hbox_pz.addWidget(w)
            hbox_pz.setAlignment(w, Qt.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ self.rb_t0, QLabel('t_x,y,z'), self.ef_tx, self.ef_ty, self.ef_tz, QLabel('q_w,x,y,z'), self.ef_qw, self.ef_qx, self.ef_qy, self.ef_qz]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_Rx)
        vbox.addLayout(hbox_Ry)
        vbox.addLayout(hbox_Rz)
        vbox.addLayout(hbox_px)
        vbox.addLayout(hbox_py)
        vbox.addLayout(hbox_pz)
        vbox.addLayout(hbox_rb)
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Transformations")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_file_action, None, quit_action))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()

if __name__ == "__main__":
    main()

