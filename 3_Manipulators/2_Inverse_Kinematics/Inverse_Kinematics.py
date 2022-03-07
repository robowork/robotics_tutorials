import sys, os, random
import numpy as np

from transforms3d import *  #supersedes deprecated transformations.py

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from mpl_toolkits import mplot3d


class RigidBodyDynamics:  

    @staticmethod
    def SkewSymmetric(p):
        p_skew = np.array([[ 0      , -p[2][0],  p[1][0]], \
                           [ p[2][0],        0, -p[0][0]], \
                           [-p[1][0],  p[0][0],        0]])
        return p_skew

    @staticmethod
    def Vee(p_skew):
        p_vee = np.array([[p_skew[2][1]], \
                          [p_skew[0][2]], \
                          [p_skew[1][0]]])
        return p_vee;

    @staticmethod
    def Inverse(T):
        R_inv = np.transpose(T[0:3,0:3])
        p = np.array([[T[0][3]], \
                      [T[1][3]], \
                      [T[2][3]]])
        T_inv = np.concatenate((np.concatenate((R_inv, -R_inv.dot(p)), axis=1), \
                                np.array([[0, 0, 0, 1]])), axis=0)
        return T_inv

    @staticmethod
    def AdjointMap(T):
        R = T[0:3,0:3]
        p = np.array([[T[0][3]], \
                      [T[1][3]], \
                      [T[2][3]]])
        p_skew = RigidBodyDynamics.SkewSymmetric(p)
        Ad_T = np.concatenate((np.concatenate((R, np.zeros((3,3))), axis=1), \
                               np.concatenate((p_skew.dot(R), R), axis=1)), axis=0)
        return Ad_T

    @staticmethod
    def ExpMap(Screw, thetaScrew_deg):
        thetaScrew = np.deg2rad(thetaScrew_deg)
        if abs(thetaScrew) < 1e-3:
            e_ScrewMatrix_thetaScrew = np.eye(4)
        elif np.linalg.norm(Screw[0:3]) < 1e-6:
            V = np.eye(3);
            e_ScrewMatrix_thetaScrew = np.concatenate((np.concatenate((np.eye(3), V.dot(Screw[3:6]) * thetaScrew), axis=1), \
                                                       np.array([[0, 0, 0, 1]])), axis=0)
        else:
            omegaSkew = RigidBodyDynamics.SkewSymmetric(Screw) * thetaScrew
            theta = np.linalg.norm(Screw[0:3] * thetaScrew)
            e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
            V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
            e_ScrewMatrix_thetaScrew = np.concatenate((np.concatenate((e_omegaSkew, V.dot(Screw[3:6] * thetaScrew)), axis=1), \
                                                       np.array([[0, 0, 0, 1]])), axis=0)
        return e_ScrewMatrix_thetaScrew

    @staticmethod
    def LogMap(T):
        R = T[0:3,0:3]
        p = np.array([[T[0][3]], \
                      [T[1][3]], \
                      [T[2][3]]])

        cos_theta = (np.trace(R)-1.0)/2
        if cos_theta > 0.999999:
            omega = np.array([[0], \
                              [0], \
                              [0]])

            V = np.eye(3);
        else:
            sin_theta = np.sqrt(1.0-cos_theta**2)
            theta = np.arccos(cos_theta)
            ln_R = (theta / (2*sin_theta)) * (R - np.transpose(R))
            omega = RigidBodyDynamics.Vee(ln_R)

            omegaSkew = ln_R
            V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)

        t = np.linalg.inv(V).dot(p)
        
        Twist = np.array([omega[0], \
                          omega[1], \
                          omega[2], \
                          t[0], \
                          t[1], \
                          t[2]])
        return Twist


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
        self.setWindowTitle('Manipulation: Inverse Kinematics')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        qApp.installEventFilter(self)

        self.fire_event = False
        self.G_T_active = np.NAN

        self.val_theta_1 = 0
        self.val_theta_2 = 0
        self.val_theta_3 = 0
        self.val_theta_4 = 0
        self.val_theta_5 = 0
        self.val_theta_6 = 0
        self.val_vel_x = 0
        self.val_vel_y = 0
        self.val_vel_z = 0
        self.val_omega_x = 0
        self.val_omega_y = 0
        self.val_omega_z = 0

        self.timer_period = 0.10
        self.timer_value = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_timer)
        self.timer.start(self.timer_period * 1000)

        self.scenario = 0

        self.cbx_theta_1.setEnabled(True)
        self.sld_theta_1.setEnabled(True)
        self.cbx_theta_2.setEnabled(True)
        self.sld_theta_2.setEnabled(True)
        self.cbx_theta_3.setEnabled(True)
        self.sld_theta_3.setEnabled(True)  
        self.cbx_theta_4.setEnabled(False)
        self.sld_theta_4.setEnabled(False)
        self.cbx_theta_5.setEnabled(False)
        self.sld_theta_5.setEnabled(False)
        self.cbx_theta_6.setEnabled(False)
        self.sld_theta_6.setEnabled(False)         
        self.cbx_vx.setEnabled(False)
        self.sld_vx.setEnabled(False)
        self.cbx_vy.setEnabled(False)
        self.sld_vy.setEnabled(False)
        self.cbx_vz.setEnabled(False)
        self.sld_vz.setEnabled(False)
        self.cbx_wx.setEnabled(False)
        self.sld_wx.setEnabled(False)
        self.cbx_wy.setEnabled(False)
        self.sld_wy.setEnabled(False)
        self.cbx_wz.setEnabled(False)
        self.sld_wz.setEnabled(False)
        self.cbx_posx.setEnabled(False)
        self.sld_posx.setEnabled(False)
        self.cbx_posy.setEnabled(False)
        self.sld_posy.setEnabled(False)
        self.cbx_posz.setEnabled(False)
        self.sld_posz.setEnabled(False)
        self.cbx_rotx.setEnabled(False)
        self.sld_rotx.setEnabled(False)
        self.cbx_roty.setEnabled(False)
        self.sld_roty.setEnabled(False)
        self.cbx_rotz.setEnabled(False)
        self.sld_rotz.setEnabled(False)

        # Space Frame / origin
        self.S_p = np.array([[0], \
                             [0], \
                             [0]])
        self.S_R = np.array([[1, 0, 0], \
                             [0, 1, 0], \
                             [0, 0, 1]])
        self.S_T = np.concatenate((np.concatenate((self.S_R, self.S_p), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
   
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

    def eventFilter(self, source, event):
        if event.type() == QEvent.KeyPress:
            #print('KeyPress: %s [%r]' % (event.key(), source))
            if event.key() == 32:  # spacebar
                self.fire_event = True
        return super(AppForm, self).eventFilter(source, event)

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
        if self.timer_value > 2:
            self.timer_value = 0
        #print(self.timer_value)

        self.on_draw()

    def on_draw(self):
        self.axes.clear()        
        self.axes.grid(True)

        if self.scenario == 0 or self.scenario == 1 or self.scenario == 2:
            # 3-Joint SE(2)
            L1 = 0.75
            L2 = 0.5
            L3 = 0.25

            Home_T = np.concatenate((np.concatenate((np.eye(3), np.array([[L1+L2+L3], [0], [0]])), axis=1), \
                                     np.array([[0, 0, 0, 1]])), axis=0)

            S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                 np.cross(-np.array([[0], [0], [1]]) , np.array([[-(L1+L2+L3)], [0], [0]]), axis=0)), axis=0)
            S2 = np.concatenate((np.array([[0], [0], [0]]), \
                                 np.array([[1], [0], [0]])), axis=0)
            S3 = np.concatenate((np.array([[0], [0], [1]]), \
                                 np.cross(-np.array([[0], [0], [1]]) , np.array([[-(L3)], [0], [0]]), axis=0)), axis=0)

            e_S1Matrix_theta1 = RigidBodyDynamics.ExpMap(S1, self.val_theta_1)
            e_S2Matrix_theta2 = RigidBodyDynamics.ExpMap(S2, self.val_theta_2/np.pi)
            e_S3Matrix_theta3 = RigidBodyDynamics.ExpMap(S3, self.val_theta_3)

            E_T = ((Home_T.dot(e_S1Matrix_theta1)).dot(e_S2Matrix_theta2)).dot(e_S3Matrix_theta3)  # Post-multiply

            # Intermediate frames
            # Joint 1
            Home_j1 = np.concatenate((np.concatenate((np.eye(3), np.array([[(0)], [0], [0]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)
             
            j1_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(0)], [0], [0]]), axis=0)), axis=0)
           
            e_j1S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j1_S1, self.val_theta_1)

            j1_T = Home_j1.dot(e_j1S1Matrix_theta1)  # Post-multiply

            # Joint 2
            Home_j2 = np.concatenate((np.concatenate((np.eye(3), np.array([[(L1+L2)], [0], [0]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j2_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(L1+L2)], [0], [0]]), axis=0)), axis=0)
            j2_S2 = np.concatenate((np.array([[0], [0], [0]]), \
                                    np.array([[1], [0], [0]])), axis=0)
           
            e_j2S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j2_S1, self.val_theta_1)
            e_j2S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j2_S2, self.val_theta_2/np.pi)
  
            j2_T = (Home_j2.dot(e_j2S1Matrix_theta1)).dot(e_j2S2Matrix_theta2)  # Post-multiply
            
            # Joint 3
            Home_j3 = np.concatenate((np.concatenate((np.eye(3), np.array([[(L1+L2)], [0], [0]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j3_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(L1+L2)], [0], [0]]), axis=0)), axis=0)
            j3_S2 = np.concatenate((np.array([[0], [0], [0]]), \
                                    np.array([[1], [0], [0]])), axis=0)
            j3_S3 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(0)], [0], [0]]), axis=0)), axis=0)     

            e_j3S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j3_S1, self.val_theta_1)
            e_j3S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j3_S2, self.val_theta_2/np.pi)
            e_j3S3Matrix_theta3 = RigidBodyDynamics.ExpMap(j3_S3, self.val_theta_3)

            j3_T = ((Home_j3.dot(e_j3S1Matrix_theta1)).dot(e_j3S2Matrix_theta2)).dot(e_j3S3Matrix_theta3)  # Post-multiply

            self.quiver_Ep = E_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Ex = E_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_Ey = E_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Ez = E_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J1p = j1_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J1x = j1_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J1y = j1_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J1z = j1_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J2p = j2_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J2x = j2_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J2y = j2_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J2z = j2_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J3p = j3_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J3x = j3_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J3y = j3_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J3z = j3_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # links
            L1_point0 = np.array([[0], [0], [0]])
            L1_point1 = np.array([[+L1], [0], [0]])
            L1_point0_transformed = j1_T.dot(np.concatenate((L1_point0, np.array([[1]])), axis=0))
            L1_point1_transformed = j1_T.dot(np.concatenate((L1_point1, np.array([[1]])), axis=0))
            L1_x = [L1_point0_transformed[0][0], L1_point1_transformed[0][0]]
            L1_y = [L1_point0_transformed[1][0], L1_point1_transformed[1][0]]
            L1_z = [L1_point0_transformed[2][0], L1_point1_transformed[2][0]]
            self.axes.plot(L1_x, L1_y, L1_z, color=(0,0,0,1.0), linewidth=3.5)

            L2_point0 = np.array([[0], [0], [0]])
            L2_point1 = np.array([[-L2 -(np.deg2rad(self.val_theta_2)/np.pi)], [0], [0]])
            L2_point0_transformed = j2_T.dot(np.concatenate((L2_point0, np.array([[1]])), axis=0))
            L2_point1_transformed = j2_T.dot(np.concatenate((L2_point1, np.array([[1]])), axis=0))
            L2_x = [L2_point0_transformed[0][0], L2_point1_transformed[0][0]]
            L2_y = [L2_point0_transformed[1][0], L2_point1_transformed[1][0]]
            L2_z = [L2_point0_transformed[2][0], L2_point1_transformed[2][0]]
            self.axes.plot(L2_x, L2_y, L2_z, color=(0.5,0.5,0.5,1.0), linewidth=3.5)

            L3_point0 = np.array([[0], [0], [0]])
            L3_point1 = np.array([[+L3], [0], [0]])
            L3_point0_transformed = j3_T.dot(np.concatenate((L3_point0, np.array([[1]])), axis=0))
            L3_point1_transformed = j3_T.dot(np.concatenate((L3_point1, np.array([[1]])), axis=0))
            L3_x = [L3_point0_transformed[0][0], L3_point1_transformed[0][0]]
            L3_y = [L3_point0_transformed[1][0], L3_point1_transformed[1][0]]
            L3_z = [L3_point0_transformed[2][0], L3_point1_transformed[2][0]]
            self.axes.plot(L3_x, L3_y, L3_z, color=(0,0,0,1.0), linewidth=3.5)

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J1p[0], self.quiver_J1p[1], self.quiver_J1p[2], scale_medium*(self.quiver_J1x[0]-self.quiver_J1p[0]), scale_medium*(self.quiver_J1x[1]-self.quiver_J1p[1]), scale_medium*(self.quiver_J1x[2]-self.quiver_J1p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J1p[0], self.quiver_J1p[1], self.quiver_J1p[2], scale_medium*(self.quiver_J1y[0]-self.quiver_J1p[0]), scale_medium*(self.quiver_J1y[1]-self.quiver_J1p[1]), scale_medium*(self.quiver_J1y[2]-self.quiver_J1p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J1p[0], self.quiver_J1p[1], self.quiver_J1p[2], scale_medium*(self.quiver_J1z[0]-self.quiver_J1p[0]), scale_medium*(self.quiver_J1z[1]-self.quiver_J1p[1]), scale_medium*(self.quiver_J1z[2]-self.quiver_J1p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J2p[0], self.quiver_J2p[1], self.quiver_J2p[2], scale_medium*(self.quiver_J2x[0]-self.quiver_J2p[0]), scale_medium*(self.quiver_J2x[1]-self.quiver_J2p[1]), scale_medium*(self.quiver_J2x[2]-self.quiver_J2p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J2p[0], self.quiver_J2p[1], self.quiver_J2p[2], scale_medium*(self.quiver_J2y[0]-self.quiver_J2p[0]), scale_medium*(self.quiver_J2y[1]-self.quiver_J2p[1]), scale_medium*(self.quiver_J2y[2]-self.quiver_J2p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J2p[0], self.quiver_J2p[1], self.quiver_J2p[2], scale_medium*(self.quiver_J2z[0]-self.quiver_J2p[0]), scale_medium*(self.quiver_J2z[1]-self.quiver_J2p[1]), scale_medium*(self.quiver_J2z[2]-self.quiver_J2p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J3p[0], self.quiver_J3p[1], self.quiver_J3p[2], scale_medium*(self.quiver_J3x[0]-self.quiver_J3p[0]), scale_medium*(self.quiver_J3x[1]-self.quiver_J3p[1]), scale_medium*(self.quiver_J3x[2]-self.quiver_J3p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J3p[0], self.quiver_J3p[1], self.quiver_J3p[2], scale_medium*(self.quiver_J3y[0]-self.quiver_J3p[0]), scale_medium*(self.quiver_J3y[1]-self.quiver_J3p[1]), scale_medium*(self.quiver_J3y[2]-self.quiver_J3p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J3p[0], self.quiver_J3p[1], self.quiver_J3p[2], scale_medium*(self.quiver_J3z[0]-self.quiver_J3p[0]), scale_medium*(self.quiver_J3z[1]-self.quiver_J3p[1]), scale_medium*(self.quiver_J3z[2]-self.quiver_J3p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], scale_full*(self.quiver_Ex[0]-self.quiver_Ep[0]), scale_full*(self.quiver_Ex[1]-self.quiver_Ep[1]), scale_full*(self.quiver_Ex[2]-self.quiver_Ep[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], scale_full*(self.quiver_Ey[0]-self.quiver_Ep[0]), scale_full*(self.quiver_Ey[1]-self.quiver_Ep[1]), scale_full*(self.quiver_Ey[2]-self.quiver_Ep[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], scale_full*(self.quiver_Ez[0]-self.quiver_Ep[0]), scale_full*(self.quiver_Ez[1]-self.quiver_Ep[1]), scale_full*(self.quiver_Ez[2]-self.quiver_Ep[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-2.0, 2.0)
            self.axes.set_ylim(-2.0, 2.0)
            self.axes.set_zlim(-2.0, 2.0)

            # Jacobian
            J_3_E = S3
            J_2_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( e_S3Matrix_theta3 )).dot(S2)
            J_1_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( e_S2Matrix_theta2.dot(e_S3Matrix_theta3) )).dot(S1)
            J_E = np.concatenate((J_1_E, J_2_E, J_3_E), axis=1)

            #if np.linalg.cond( J_E.dot(np.transpose(J_E)) ) < 1/sys.float_info.epsilon:
            #    print("Non Singular")
            #else:
            #    print("Singular")

            J_E_angular = J_E[0:3,:] 
            J_E_linear = J_E[3:6,:]

            A_manip_angular = J_E_angular.dot(np.transpose(J_E_angular))
            A_rank_angular = np.linalg.matrix_rank(A_manip_angular)
            A_manip_linear = J_E_linear.dot(np.transpose(J_E_linear))
            A_rank_linear = np.linalg.matrix_rank(A_manip_linear)
            print("A_rank_linear = " + str(A_rank_linear) + " , A_rank_angular = " + str(A_rank_angular))

            W_linear, V_linear = np.linalg.eig(A_manip_linear)
            idx_linear = W_linear.argsort()[::-1]   
            W = W_linear[idx_linear]
            V = V_linear[:,idx_linear]
            #for e in range(0,len(W_linear)):
            for e in range(0,A_rank_linear):
                eVec = V[:,e]
                eVal = W[e]
                axis_halflen = np.sqrt(eVal)
   
                self.quiver_Ee = np.array([[eVec[0]], \
                                           [eVec[1]], \
                                           [eVec[2]]])
                self.quiver_Ee = E_T.dot(np.concatenate((self.quiver_Ee, np.array([[1]])), axis=0))

                self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], axis_halflen*(self.quiver_Ee[0]-self.quiver_Ep[0]), axis_halflen*(self.quiver_Ee[1]-self.quiver_Ep[1]), axis_halflen*(self.quiver_Ee[2]-self.quiver_Ep[2]), color=['m'], arrow_length_ratio=0.15)
                self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], -axis_halflen*(self.quiver_Ee[0]-self.quiver_Ep[0]), -axis_halflen*(self.quiver_Ee[1]-self.quiver_Ep[1]), -axis_halflen*(self.quiver_Ee[2]-self.quiver_Ep[2]), color=['m'], arrow_length_ratio=0.15)

            ellipsoid_u = np.linspace(0, 2*np.pi, 25)
            ellipsoid_v = np.linspace(0, np.pi, 25)
            ellipsoid_x = np.sqrt(W[0]) * np.outer(np.cos(ellipsoid_u), np.sin(ellipsoid_v))
            ellipsoid_y = np.sqrt(W[1]) * np.outer(np.sin(ellipsoid_u), np.sin(ellipsoid_v))
            ellipsoid_z = np.sqrt(W[2]) * np.outer(np.ones_like(ellipsoid_u), np.cos(ellipsoid_v))
            ellipsoid_x_transformed = ellipsoid_x.copy()
            ellipsoid_y_transformed = ellipsoid_y.copy()
            ellipsoid_z_transformed = ellipsoid_z.copy()
            for j in range(len(ellipsoid_x_transformed)):
                for i in range(len(ellipsoid_x_transformed[j])):
                    ellipsoid_i_xyz = E_T.dot(np.concatenate((np.concatenate((V,np.array([[0],[0],[0]])), axis=1), np.array([[0,0,0,1]])), axis=0)).dot(np.array([[ellipsoid_x_transformed[i][j]], [ellipsoid_y_transformed[i][j]], [ellipsoid_z_transformed[i][j]], [1]]))
                    ellipsoid_x_transformed[i][j] = ellipsoid_i_xyz[0][0]
                    ellipsoid_y_transformed[i][j] = ellipsoid_i_xyz[1][0]
                    ellipsoid_z_transformed[i][j] = ellipsoid_i_xyz[2][0]
            self.axes.plot_surface(ellipsoid_x_transformed, ellipsoid_y_transformed, ellipsoid_z_transformed, rstride=4, cstride=4, cmap='seismic', alpha=0.125)

        if self.scenario == 1:
            # Use Jacobian
            J_S = RigidBodyDynamics.AdjointMap(E_T).dot(J_E)

            J_S_pinv = np.linalg.pinv(J_S)
            
            Twist_S = np.array([[np.deg2rad(self.val_omega_x)], \
                                [np.deg2rad(self.val_omega_y)], \
                                [np.deg2rad(self.val_omega_z)], \
                                [self.val_vel_x], \
                                [self.val_vel_y], \
                                [self.val_vel_z]])     
       
            Joints_vel = J_S_pinv.dot(Twist_S)

            # revolute
            self.val_theta_1 = self.val_theta_1 + np.rad2deg(Joints_vel[0]) * 0.1
            if self.val_theta_1 > 180.0:
                self.val_theta_1 = 180.0
            elif self.val_theta_1 < -180.0:
                self.val_theta_1 = -180.0

            # prismatic
            self.val_theta_2 = self.val_theta_2 + np.rad2deg(Joints_vel[1]) * 0.1
            if self.val_theta_2 > 180.0:
                self.val_theta_2 = 180.0
            elif self.val_theta_2 < -180.0:
                self.val_theta_2 = -180.0

            # wrap: revolute
            self.val_theta_3 = self.val_theta_3 + np.rad2deg(Joints_vel[2]) * 0.1
            if self.val_theta_3 > 180.0:
                self.val_theta_3 = 180.0
            elif self.val_theta_3 < -180.0:
                self.val_theta_3 = -180.0

            self.sld_theta_1.blockSignals(True)
            self.sld_theta_1.setValue(self.val_theta_1)
            self.sld_theta_1.blockSignals(False)

            self.sld_theta_2.blockSignals(True)
            self.sld_theta_2.setValue(self.val_theta_2)
            self.sld_theta_2.blockSignals(False)

            self.sld_theta_3.blockSignals(True)
            self.sld_theta_3.setValue(self.val_theta_3)
            self.sld_theta_3.blockSignals(False)

        if self.scenario == 2:
            # Goal transformation (Body Frame-based)
            G_p_zyx = np.array([[self.val_pos_x], \
                                [self.val_pos_y], \
                                [self.val_pos_z]])
            G_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), G_p_zyx), axis=1), \
                                            np.array([[0, 0, 0, 1]])), axis=0)
            angle_z = np.deg2rad(self.val_rot_z)
            angle_y = np.deg2rad(self.val_rot_y)
            angle_x = np.deg2rad(self.val_rot_x)
            G_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0,                0,               1]])
            G_T_z = np.concatenate((np.concatenate((G_R_z, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            G_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            G_T_y = np.concatenate((np.concatenate((G_R_y, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            G_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
            G_T_x = np.concatenate((np.concatenate((G_R_x, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)

            G_T_zyx = (((self.S_T.dot(G_T_zyx_trans)).dot(G_T_z)).dot(G_T_y)).dot(G_T_x)

            G_T = G_T_zyx
            G_p = G_T_zyx[0:3,3]
            G_R = G_T_zyx[0:3,0:3]

            self.quiver_Gp = G_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Gx = G_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_Gy = G_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Gz = G_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_goal = 1.5
            self.axes.quiver(self.quiver_Gp[0], self.quiver_Gp[1], self.quiver_Gp[2], scale_goal*(self.quiver_Gx[0]-self.quiver_Gp[0]), scale_goal*(self.quiver_Gx[1]-self.quiver_Gp[1]), scale_goal*(self.quiver_Gx[2]-self.quiver_Gp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Gp[0], self.quiver_Gp[1], self.quiver_Gp[2], scale_goal*(self.quiver_Gy[0]-self.quiver_Gp[0]), scale_goal*(self.quiver_Gy[1]-self.quiver_Gp[1]), scale_goal*(self.quiver_Gy[2]-self.quiver_Gp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Gp[0], self.quiver_Gp[1], self.quiver_Gp[2], scale_goal*(self.quiver_Gz[0]-self.quiver_Gp[0]), scale_goal*(self.quiver_Gz[1]-self.quiver_Gp[1]), scale_goal*(self.quiver_Gz[2]-self.quiver_Gp[2]), color=['b'], arrow_length_ratio=0.15)

            # update goal for IK problem
            if self.fire_event:
                self.fire_event = False
 
                self.G_T_active = G_T

            if not np.isnan(self.G_T_active).any():
                # Solve IK problem at i-th iteration
                error_i = RigidBodyDynamics.Inverse(E_T).dot(self.G_T_active) 
                Twist_i = RigidBodyDynamics.LogMap(error_i)

                # Use Jacobian
                J_E_pinv = np.linalg.pinv(J_E) 
       
                Joints_vel = J_E_pinv.dot(Twist_i)

                # revolute
                self.val_theta_1 = self.val_theta_1 + np.rad2deg(Joints_vel[0]) * 0.1
                if self.val_theta_1 > 180.0:
                    self.val_theta_1 = 180.0
                elif self.val_theta_1 < -180.0:
                    self.val_theta_1 = -180.0

                # revolute
                self.val_theta_2 = self.val_theta_2 + np.rad2deg(Joints_vel[1]) * 0.1
                if self.val_theta_2 > 180.0:
                    self.val_theta_2 = 180.0
                elif self.val_theta_2 < -180.0:
                    self.val_theta_2 = -180.0

                # revolute
                self.val_theta_3 = self.val_theta_3 + np.rad2deg(Joints_vel[2]) * 0.1
                if self.val_theta_3 > 180.0:
                    self.val_theta_3 = 180.0
                elif self.val_theta_3 < -180.0:
                    self.val_theta_3 = -180.0

                self.sld_theta_1.blockSignals(True)
                self.sld_theta_1.setValue(self.val_theta_1)
                self.sld_theta_1.blockSignals(False)

                self.sld_theta_2.blockSignals(True)
                self.sld_theta_2.setValue(self.val_theta_2)
                self.sld_theta_2.blockSignals(False)

                self.sld_theta_3.blockSignals(True)
                self.sld_theta_3.setValue(self.val_theta_3)
                self.sld_theta_3.blockSignals(False)

        if self.scenario == 3 or self.scenario == 4 or self.scenario == 5:
            # 6-Joint SE(3)
            L1 = 0.50
            L2 = 0.25
            L3 = 1.0
            L4 = 0.50
            L5 = 1.0
            L6 = 0.25
            L7 = 0.50
            L8 = 0.50
            L9 = 0.25
            L10 = 0.10

            Home_T = np.concatenate((np.concatenate((np.eye(3), np.array([[L9+L10], [-L2+L4-L6], [L1+L3+L5+L7+L8]])), axis=1), \
                                     np.array([[0, 0, 0, 1]])), axis=0)

            S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                 np.cross(-np.array([[0], [0], [1]]) , np.array([[-L10-L9], [+L6-L4+L2], [-L8-L7-L5-L3-L1]]), axis=0)), axis=0)
            S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                 np.cross(-np.array([[0], [1], [0]]) , np.array([[-L10-L9], [+L6-L4], [-L8-L7-L5-L3]]), axis=0)), axis=0)
            S3 = np.concatenate((np.array([[0], [1], [0]]), \
                                 np.cross(-np.array([[0], [1], [0]]) , np.array([[-L10-L9], [+L6], [-L8-L7-L5]]), axis=0)), axis=0)
            S4 = np.concatenate((np.array([[0], [1], [0]]), \
                                 np.cross(-np.array([[0], [1], [0]]) , np.array([[-L10-L9], [0], [-L8-L7]]), axis=0)), axis=0)
            S5 = np.concatenate((np.array([[0], [0], [1]]), \
                                 np.cross(-np.array([[0], [0], [1]]) , np.array([[-L10-L9], [0], [-L8]]), axis=0)), axis=0)
            S6 = np.concatenate((np.array([[1], [0], [0]]), \
                                 np.cross(-np.array([[1], [0], [0]]) , np.array([[-L10], [0], [0]]), axis=0)), axis=0)

            e_S1Matrix_theta1 = RigidBodyDynamics.ExpMap(S1, self.val_theta_1)
            e_S2Matrix_theta2 = RigidBodyDynamics.ExpMap(S2, self.val_theta_2)
            e_S3Matrix_theta3 = RigidBodyDynamics.ExpMap(S3, self.val_theta_3)
            e_S4Matrix_theta4 = RigidBodyDynamics.ExpMap(S4, self.val_theta_4)
            e_S5Matrix_theta5 = RigidBodyDynamics.ExpMap(S5, self.val_theta_5)
            e_S6Matrix_theta6 = RigidBodyDynamics.ExpMap(S6, self.val_theta_6)

            E_T = (((((Home_T.dot(e_S1Matrix_theta1)).dot(e_S2Matrix_theta2)).dot(e_S3Matrix_theta3)).dot(e_S4Matrix_theta4)).dot(e_S5Matrix_theta5)).dot(e_S6Matrix_theta6)  # Post-multiply

            # Intermediate frames
            # Joint 1
            Home_j1 = np.concatenate((np.concatenate((np.eye(3), np.array([[0], [0], [0]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)
             
            j1_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[0], [0], [0]]), axis=0)), axis=0)
           
            e_j1S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j1_S1, self.val_theta_1)

            j1_T = Home_j1.dot(e_j1S1Matrix_theta1)  # Post-multiply

            # Joint 2
            Home_j2 = np.concatenate((np.concatenate((np.eye(3), np.array([[0], [-L2], [L1]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j2_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[0], [+L2], [-L1]]), axis=0)), axis=0)
            j2_S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [0], [0]]), axis=0)), axis=0)
           
            e_j2S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j2_S1, self.val_theta_1)
            e_j2S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j2_S2, self.val_theta_2)

            j2_T = (Home_j2.dot(e_j2S1Matrix_theta1)).dot(e_j2S2Matrix_theta2)  # Post-multiply
            
            # Joint 3
            Home_j3 = np.concatenate((np.concatenate((np.eye(3), np.array([[0], [-L2+L4], [L1+L3]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j3_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[0], [-L4+L2], [-L3-L1]]), axis=0)), axis=0)
            j3_S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [-L4], [-L3]]), axis=0)), axis=0)   
            j3_S3 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [0], [0]]), axis=0)), axis=0) 

            e_j3S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j3_S1, self.val_theta_1)
            e_j3S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j3_S2, self.val_theta_2)
            e_j3S3Matrix_theta3 = RigidBodyDynamics.ExpMap(j3_S3, self.val_theta_3)

            j3_T = ((Home_j3.dot(e_j3S1Matrix_theta1)).dot(e_j3S2Matrix_theta2)).dot(e_j3S3Matrix_theta3)  # Post-multiply

            # Joint 4
            Home_j4 = np.concatenate((np.concatenate((np.eye(3), np.array([[0], [-L2+L4-L6], [L1+L3+L5]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j4_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[0], [+L6-L4+L2], [-L5-L3-L1]]), axis=0)), axis=0)
            j4_S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [+L6-L4], [-L5-L3]]), axis=0)), axis=0)  
            j4_S3 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [+L6], [-L5]]), axis=0)), axis=0) 
            j4_S4 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [0], [0]]), axis=0)), axis=0) 

            e_j4S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j4_S1, self.val_theta_1)
            e_j4S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j4_S2, self.val_theta_2)
            e_j4S3Matrix_theta3 = RigidBodyDynamics.ExpMap(j4_S3, self.val_theta_3)
            e_j4S4Matrix_theta4 = RigidBodyDynamics.ExpMap(j4_S4, self.val_theta_4)

            j4_T = (((Home_j4.dot(e_j4S1Matrix_theta1)).dot(e_j4S2Matrix_theta2)).dot(e_j4S3Matrix_theta3)).dot(e_j4S4Matrix_theta4)  # Post-multiply

            # Joint 5
            Home_j5 = np.concatenate((np.concatenate((np.eye(3), np.array([[0], [-L2+L4-L6], [L1+L3+L5+L7]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j5_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[0], [+L6-L4+L2], [-L7-L5-L3-L1]]), axis=0)), axis=0)
            j5_S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [+L6-L4], [-L7-L5-L3]]), axis=0)), axis=0)  
            j5_S3 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [+L6], [-L7-L5]]), axis=0)), axis=0) 
            j5_S4 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [0], [-L7]]), axis=0)), axis=0) 
            j5_S5 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[0], [0], [0]]), axis=0)), axis=0) 

            e_j5S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j5_S1, self.val_theta_1)
            e_j5S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j5_S2, self.val_theta_2)
            e_j5S3Matrix_theta3 = RigidBodyDynamics.ExpMap(j5_S3, self.val_theta_3)
            e_j5S4Matrix_theta4 = RigidBodyDynamics.ExpMap(j5_S4, self.val_theta_4)
            e_j5S5Matrix_theta5 = RigidBodyDynamics.ExpMap(j5_S5, self.val_theta_5)

            j5_T = ((((Home_j5.dot(e_j5S1Matrix_theta1)).dot(e_j5S2Matrix_theta2)).dot(e_j5S3Matrix_theta3)).dot(e_j5S4Matrix_theta4)).dot(e_j5S5Matrix_theta5)  # Post-multiply

            # Joint 6
            Home_j6 = np.concatenate((np.concatenate((np.eye(3), np.array([[L9], [-L2+L4-L6], [L1+L3+L5+L7+L8]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j6_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-L9], [+L6-L4+L2], [-L8-L7-L5-L3-L1]]), axis=0)), axis=0)
            j6_S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[-L9], [+L6-L4], [-L8-L7-L5-L3]]), axis=0)), axis=0)  
            j6_S3 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[-L9], [+L6], [-L8-L7-L5]]), axis=0)), axis=0) 
            j6_S4 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[-L9], [0], [-L8-L7]]), axis=0)), axis=0) 
            j6_S5 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-L9], [0], [-L8]]), axis=0)), axis=0) 
            j6_S6 = np.concatenate((np.array([[1], [0], [0]]), \
                                    np.cross(-np.array([[1], [0], [0]]) , np.array([[0], [0], [0]]), axis=0)), axis=0) 

            e_j6S1Matrix_theta1 = RigidBodyDynamics.ExpMap(j6_S1, self.val_theta_1)
            e_j6S2Matrix_theta2 = RigidBodyDynamics.ExpMap(j6_S2, self.val_theta_2)
            e_j6S3Matrix_theta3 = RigidBodyDynamics.ExpMap(j6_S3, self.val_theta_3)
            e_j6S4Matrix_theta4 = RigidBodyDynamics.ExpMap(j6_S4, self.val_theta_4)
            e_j6S5Matrix_theta5 = RigidBodyDynamics.ExpMap(j6_S5, self.val_theta_5)
            e_j6S6Matrix_theta6 = RigidBodyDynamics.ExpMap(j6_S6, self.val_theta_6)

            j6_T = (((((Home_j6.dot(e_j6S1Matrix_theta1)).dot(e_j6S2Matrix_theta2)).dot(e_j6S3Matrix_theta3)).dot(e_j6S4Matrix_theta4)).dot(e_j6S5Matrix_theta5)).dot(e_j6S6Matrix_theta6)  # Post-multiply

            self.quiver_Ep = E_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Ex = E_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_Ey = E_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Ez = E_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J1p = j1_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J1x = j1_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J1y = j1_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J1z = j1_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J2p = j2_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J2x = j2_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J2y = j2_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J2z = j2_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J3p = j3_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J3x = j3_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J3y = j3_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J3z = j3_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J4p = j4_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J4x = j4_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J4y = j4_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J4z = j4_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J5p = j5_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J5x = j5_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J5y = j5_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J5z = j5_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            self.quiver_J6p = j6_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_J6x = j6_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_J6y = j6_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_J6z = j6_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # links
            L1_point0 = np.array([[0], [0], [0]])
            L1_point1 = np.array([[0], [0], [+L1]])
            L1_point2 = np.array([[0], [-L2], [+L1]])
            L1_point0_transformed = j1_T.dot(np.concatenate((L1_point0, np.array([[1]])), axis=0))
            L1_point1_transformed = j1_T.dot(np.concatenate((L1_point1, np.array([[1]])), axis=0))
            L1_point2_transformed = j1_T.dot(np.concatenate((L1_point2, np.array([[1]])), axis=0))
            L1_x = [L1_point0_transformed[0][0], L1_point1_transformed[0][0], L1_point2_transformed[0][0]]
            L1_y = [L1_point0_transformed[1][0], L1_point1_transformed[1][0], L1_point2_transformed[1][0]]
            L1_z = [L1_point0_transformed[2][0], L1_point1_transformed[2][0], L1_point2_transformed[2][0]]
            self.axes.plot(L1_x, L1_y, L1_z, color=(0,0,0,1.0), linewidth=3.5)

            L2_point0 = np.array([[0], [0], [0]])
            L2_point1 = np.array([[0], [0], [+L3]])
            L2_point2 = np.array([[0], [+L4], [+L3]])
            L2_point0_transformed = j2_T.dot(np.concatenate((L2_point0, np.array([[1]])), axis=0))
            L2_point1_transformed = j2_T.dot(np.concatenate((L2_point1, np.array([[1]])), axis=0))
            L2_point2_transformed = j2_T.dot(np.concatenate((L2_point2, np.array([[1]])), axis=0))
            L2_x = [L2_point0_transformed[0][0], L2_point1_transformed[0][0], L2_point2_transformed[0][0]]
            L2_y = [L2_point0_transformed[1][0], L2_point1_transformed[1][0], L2_point2_transformed[1][0]]
            L2_z = [L2_point0_transformed[2][0], L2_point1_transformed[2][0], L2_point2_transformed[2][0]]
            self.axes.plot(L2_x, L2_y, L2_z, color=(0.5,0.5,0.5,1.0), linewidth=3.5)

            L3_point0 = np.array([[0], [0], [0]])
            L3_point1 = np.array([[0], [0], [+L5]])
            L3_point2 = np.array([[0], [-L6], [+L5]])
            L3_point0_transformed = j3_T.dot(np.concatenate((L3_point0, np.array([[1]])), axis=0))
            L3_point1_transformed = j3_T.dot(np.concatenate((L3_point1, np.array([[1]])), axis=0))
            L3_point2_transformed = j3_T.dot(np.concatenate((L3_point2, np.array([[1]])), axis=0))
            L3_x = [L3_point0_transformed[0][0], L3_point1_transformed[0][0], L3_point2_transformed[0][0]]
            L3_y = [L3_point0_transformed[1][0], L3_point1_transformed[1][0], L3_point2_transformed[1][0]]
            L3_z = [L3_point0_transformed[2][0], L3_point1_transformed[2][0], L3_point2_transformed[2][0]]
            self.axes.plot(L3_x, L3_y, L3_z, color=(0,0,0,1.0), linewidth=3.5)

            L4_point0 = np.array([[0], [0], [0]])
            L4_point1 = np.array([[0], [0], [+L7]])
            L4_point0_transformed = j4_T.dot(np.concatenate((L4_point0, np.array([[1]])), axis=0))
            L4_point1_transformed = j4_T.dot(np.concatenate((L4_point1, np.array([[1]])), axis=0))
            L4_x = [L4_point0_transformed[0][0], L4_point1_transformed[0][0]]
            L4_y = [L4_point0_transformed[1][0], L4_point1_transformed[1][0]]
            L4_z = [L4_point0_transformed[2][0], L4_point1_transformed[2][0]]
            self.axes.plot(L4_x, L4_y, L4_z, color=(0.5,0.5,0.5,1.0), linewidth=3.5)

            L5_point0 = np.array([[0], [0], [0]])
            L5_point1 = np.array([[0], [0], [+L8]])
            L5_point2 = np.array([[+L9], [0], [+L8]])
            L5_point0_transformed = j5_T.dot(np.concatenate((L5_point0, np.array([[1]])), axis=0))
            L5_point1_transformed = j5_T.dot(np.concatenate((L5_point1, np.array([[1]])), axis=0))
            L5_point2_transformed = j5_T.dot(np.concatenate((L5_point2, np.array([[1]])), axis=0))
            L5_x = [L5_point0_transformed[0][0], L5_point1_transformed[0][0], L5_point2_transformed[0][0]]
            L5_y = [L5_point0_transformed[1][0], L5_point1_transformed[1][0], L5_point2_transformed[1][0]]
            L5_z = [L5_point0_transformed[2][0], L5_point1_transformed[2][0], L5_point2_transformed[2][0]]
            self.axes.plot(L5_x, L5_y, L5_z, color=(0,0,0,1.0), linewidth=3.5)

            L6_point0 = np.array([[0], [0], [0]])
            L6_point1 = np.array([[+L10], [0], [0]])
            L6_point0_transformed = j6_T.dot(np.concatenate((L6_point0, np.array([[1]])), axis=0))
            L6_point1_transformed = j6_T.dot(np.concatenate((L6_point1, np.array([[1]])), axis=0))
            L6_x = [L6_point0_transformed[0][0], L6_point1_transformed[0][0]]
            L6_y = [L6_point0_transformed[1][0], L6_point1_transformed[1][0]]
            L6_z = [L6_point0_transformed[2][0], L6_point1_transformed[2][0]]
            self.axes.plot(L6_x, L6_y, L6_z, color=(0.5,0.5,0.5,1.0), linewidth=3.5)

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J1p[0], self.quiver_J1p[1], self.quiver_J1p[2], scale_medium*(self.quiver_J1x[0]-self.quiver_J1p[0]), scale_medium*(self.quiver_J1x[1]-self.quiver_J1p[1]), scale_medium*(self.quiver_J1x[2]-self.quiver_J1p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J1p[0], self.quiver_J1p[1], self.quiver_J1p[2], scale_medium*(self.quiver_J1y[0]-self.quiver_J1p[0]), scale_medium*(self.quiver_J1y[1]-self.quiver_J1p[1]), scale_medium*(self.quiver_J1y[2]-self.quiver_J1p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J1p[0], self.quiver_J1p[1], self.quiver_J1p[2], scale_medium*(self.quiver_J1z[0]-self.quiver_J1p[0]), scale_medium*(self.quiver_J1z[1]-self.quiver_J1p[1]), scale_medium*(self.quiver_J1z[2]-self.quiver_J1p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J2p[0], self.quiver_J2p[1], self.quiver_J2p[2], scale_medium*(self.quiver_J2x[0]-self.quiver_J2p[0]), scale_medium*(self.quiver_J2x[1]-self.quiver_J2p[1]), scale_medium*(self.quiver_J2x[2]-self.quiver_J2p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J2p[0], self.quiver_J2p[1], self.quiver_J2p[2], scale_medium*(self.quiver_J2y[0]-self.quiver_J2p[0]), scale_medium*(self.quiver_J2y[1]-self.quiver_J2p[1]), scale_medium*(self.quiver_J2y[2]-self.quiver_J2p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J2p[0], self.quiver_J2p[1], self.quiver_J2p[2], scale_medium*(self.quiver_J2z[0]-self.quiver_J2p[0]), scale_medium*(self.quiver_J2z[1]-self.quiver_J2p[1]), scale_medium*(self.quiver_J2z[2]-self.quiver_J2p[2]), color=['b'], arrow_length_ratio=0.15)         
            self.axes.quiver(self.quiver_J3p[0], self.quiver_J3p[1], self.quiver_J3p[2], scale_medium*(self.quiver_J3x[0]-self.quiver_J3p[0]), scale_medium*(self.quiver_J3x[1]-self.quiver_J3p[1]), scale_medium*(self.quiver_J3x[2]-self.quiver_J3p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J3p[0], self.quiver_J3p[1], self.quiver_J3p[2], scale_medium*(self.quiver_J3y[0]-self.quiver_J3p[0]), scale_medium*(self.quiver_J3y[1]-self.quiver_J3p[1]), scale_medium*(self.quiver_J3y[2]-self.quiver_J3p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J3p[0], self.quiver_J3p[1], self.quiver_J3p[2], scale_medium*(self.quiver_J3z[0]-self.quiver_J3p[0]), scale_medium*(self.quiver_J3z[1]-self.quiver_J3p[1]), scale_medium*(self.quiver_J3z[2]-self.quiver_J3p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J4p[0], self.quiver_J4p[1], self.quiver_J4p[2], scale_medium*(self.quiver_J4x[0]-self.quiver_J4p[0]), scale_medium*(self.quiver_J4x[1]-self.quiver_J4p[1]), scale_medium*(self.quiver_J4x[2]-self.quiver_J4p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J4p[0], self.quiver_J4p[1], self.quiver_J4p[2], scale_medium*(self.quiver_J4y[0]-self.quiver_J4p[0]), scale_medium*(self.quiver_J4y[1]-self.quiver_J4p[1]), scale_medium*(self.quiver_J4y[2]-self.quiver_J4p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J4p[0], self.quiver_J4p[1], self.quiver_J4p[2], scale_medium*(self.quiver_J4z[0]-self.quiver_J4p[0]), scale_medium*(self.quiver_J4z[1]-self.quiver_J4p[1]), scale_medium*(self.quiver_J4z[2]-self.quiver_J4p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J5p[0], self.quiver_J5p[1], self.quiver_J5p[2], scale_medium*(self.quiver_J5x[0]-self.quiver_J5p[0]), scale_medium*(self.quiver_J5x[1]-self.quiver_J5p[1]), scale_medium*(self.quiver_J5x[2]-self.quiver_J5p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J5p[0], self.quiver_J5p[1], self.quiver_J5p[2], scale_medium*(self.quiver_J5y[0]-self.quiver_J5p[0]), scale_medium*(self.quiver_J5y[1]-self.quiver_J5p[1]), scale_medium*(self.quiver_J5y[2]-self.quiver_J5p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J5p[0], self.quiver_J5p[1], self.quiver_J5p[2], scale_medium*(self.quiver_J5z[0]-self.quiver_J5p[0]), scale_medium*(self.quiver_J5z[1]-self.quiver_J5p[1]), scale_medium*(self.quiver_J5z[2]-self.quiver_J5p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J6p[0], self.quiver_J6p[1], self.quiver_J6p[2], scale_medium*(self.quiver_J6x[0]-self.quiver_J6p[0]), scale_medium*(self.quiver_J6x[1]-self.quiver_J6p[1]), scale_medium*(self.quiver_J6x[2]-self.quiver_J6p[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J6p[0], self.quiver_J6p[1], self.quiver_J6p[2], scale_medium*(self.quiver_J6y[0]-self.quiver_J6p[0]), scale_medium*(self.quiver_J6y[1]-self.quiver_J6p[1]), scale_medium*(self.quiver_J6y[2]-self.quiver_J6p[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_J6p[0], self.quiver_J6p[1], self.quiver_J6p[2], scale_medium*(self.quiver_J6z[0]-self.quiver_J6p[0]), scale_medium*(self.quiver_J6z[1]-self.quiver_J6p[1]), scale_medium*(self.quiver_J6z[2]-self.quiver_J6p[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], scale_full*(self.quiver_Ex[0]-self.quiver_Ep[0]), scale_full*(self.quiver_Ex[1]-self.quiver_Ep[1]), scale_full*(self.quiver_Ex[2]-self.quiver_Ep[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], scale_full*(self.quiver_Ey[0]-self.quiver_Ep[0]), scale_full*(self.quiver_Ey[1]-self.quiver_Ep[1]), scale_full*(self.quiver_Ey[2]-self.quiver_Ep[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], scale_full*(self.quiver_Ez[0]-self.quiver_Ep[0]), scale_full*(self.quiver_Ez[1]-self.quiver_Ep[1]), scale_full*(self.quiver_Ez[2]-self.quiver_Ep[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-3.0, 3.0)
            self.axes.set_ylim(-3.0, 3.0)
            self.axes.set_zlim(-3.0, 3.0)

            # Jacobian
            J_6_E = S6
            J_5_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( e_S6Matrix_theta6 )).dot(S5)
            J_4_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( e_S5Matrix_theta5.dot(e_S6Matrix_theta6) )).dot(S4)
            J_3_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( (e_S4Matrix_theta4.dot(e_S5Matrix_theta5)).dot(e_S6Matrix_theta6) )).dot(S3)
            J_2_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( ((e_S3Matrix_theta3.dot(e_S4Matrix_theta4)).dot(e_S5Matrix_theta5)).dot(e_S6Matrix_theta6) )).dot(S2)
            J_1_E = RigidBodyDynamics.AdjointMap(RigidBodyDynamics.Inverse( (((e_S2Matrix_theta2.dot(e_S3Matrix_theta3)).dot(e_S4Matrix_theta4)).dot(e_S5Matrix_theta5)).dot(e_S6Matrix_theta6) )).dot(S1)
            J_E = np.concatenate((J_1_E, J_2_E, J_3_E, J_4_E, J_5_E, J_6_E), axis=1)

            #if np.linalg.cond( J_E.dot(np.transpose(J_E)) ) < 1/sys.float_info.epsilon:
            #    print("Non Singular")
            #else:
            #    print("Singular")

            J_E_angular = J_E[0:3,:] 
            J_E_linear = J_E[3:6,:]

            A_manip_angular = J_E_angular.dot(np.transpose(J_E_angular))
            A_rank_angular = np.linalg.matrix_rank(A_manip_angular)
            A_manip_linear = J_E_linear.dot(np.transpose(J_E_linear))
            A_rank_linear = np.linalg.matrix_rank(A_manip_linear)
            print("A_rank_linear = " + str(A_rank_linear) + " , A_rank_angular = " + str(A_rank_angular))

            W_linear, V_linear = np.linalg.eig(A_manip_linear)
            idx_linear = W_linear.argsort()[::-1]   
            W = W_linear[idx_linear]
            V = V_linear[:,idx_linear]
            #for e in range(0,len(W_linear)):
            for e in range(0,A_rank_linear):
                eVec = V[:,e]
                eVal = W[e]
                axis_halflen = np.sqrt(eVal)
   
                self.quiver_Ee = np.array([[eVec[0]], \
                                           [eVec[1]], \
                                           [eVec[2]]])
                self.quiver_Ee = E_T.dot(np.concatenate((self.quiver_Ee, np.array([[1]])), axis=0))

                self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], axis_halflen*(self.quiver_Ee[0]-self.quiver_Ep[0]), axis_halflen*(self.quiver_Ee[1]-self.quiver_Ep[1]), axis_halflen*(self.quiver_Ee[2]-self.quiver_Ep[2]), color=['m'], arrow_length_ratio=0.15)
                self.axes.quiver(self.quiver_Ep[0], self.quiver_Ep[1], self.quiver_Ep[2], -axis_halflen*(self.quiver_Ee[0]-self.quiver_Ep[0]), -axis_halflen*(self.quiver_Ee[1]-self.quiver_Ep[1]), -axis_halflen*(self.quiver_Ee[2]-self.quiver_Ep[2]), color=['m'], arrow_length_ratio=0.15)

            ellipsoid_u = np.linspace(0, 2*np.pi, 25)
            ellipsoid_v = np.linspace(0, np.pi, 25)
            ellipsoid_x = np.sqrt(W[0]) * np.outer(np.cos(ellipsoid_u), np.sin(ellipsoid_v))
            ellipsoid_y = np.sqrt(W[1]) * np.outer(np.sin(ellipsoid_u), np.sin(ellipsoid_v))
            ellipsoid_z = np.sqrt(W[2]) * np.outer(np.ones_like(ellipsoid_u), np.cos(ellipsoid_v))
            ellipsoid_x_transformed = ellipsoid_x.copy()
            ellipsoid_y_transformed = ellipsoid_y.copy()
            ellipsoid_z_transformed = ellipsoid_z.copy()
            for j in range(len(ellipsoid_x_transformed)):
                for i in range(len(ellipsoid_x_transformed[j])):
                    ellipsoid_i_xyz = E_T.dot(np.concatenate((np.concatenate((V,np.array([[0],[0],[0]])), axis=1), np.array([[0,0,0,1]])), axis=0)).dot(np.array([[ellipsoid_x_transformed[i][j]], [ellipsoid_y_transformed[i][j]], [ellipsoid_z_transformed[i][j]], [1]]))
                    ellipsoid_x_transformed[i][j] = ellipsoid_i_xyz[0][0]
                    ellipsoid_y_transformed[i][j] = ellipsoid_i_xyz[1][0]
                    ellipsoid_z_transformed[i][j] = ellipsoid_i_xyz[2][0]
            self.axes.plot_surface(ellipsoid_x_transformed, ellipsoid_y_transformed, ellipsoid_z_transformed, rstride=4, cstride=4, cmap='seismic', alpha=0.125)

        if self.scenario == 4:
            # Use Jacobian
            J_S = RigidBodyDynamics.AdjointMap(E_T).dot(J_E)

            J_S_pinv = np.linalg.pinv(J_S)
            
            Twist_S = np.array([[np.deg2rad(self.val_omega_x)], \
                                [np.deg2rad(self.val_omega_y)], \
                                [np.deg2rad(self.val_omega_z)], \
                                [self.val_vel_x], \
                                [self.val_vel_y], \
                                [self.val_vel_z]])     
       
            Joints_vel = J_S_pinv.dot(Twist_S)

            # revolute
            self.val_theta_1 = self.val_theta_1 + np.rad2deg(Joints_vel[0]) * 0.1
            if self.val_theta_1 > 180.0:
                self.val_theta_1 = 180.0
            elif self.val_theta_1 < -180.0:
                self.val_theta_1 = -180.0

            # revolute
            self.val_theta_2 = self.val_theta_2 + np.rad2deg(Joints_vel[1]) * 0.1
            if self.val_theta_2 > 180.0:
                self.val_theta_2 = 180.0
            elif self.val_theta_2 < -180.0:
                self.val_theta_2 = -180.0

            # revolute
            self.val_theta_3 = self.val_theta_3 + np.rad2deg(Joints_vel[2]) * 0.1
            if self.val_theta_3 > 180.0:
                self.val_theta_3 = 180.0
            elif self.val_theta_3 < -180.0:
                self.val_theta_3 = -180.0

            # revolute
            self.val_theta_4 = self.val_theta_4 + np.rad2deg(Joints_vel[3]) * 0.1
            if self.val_theta_4 > 180.0:
                self.val_theta_4 = 180.0
            elif self.val_theta_4 < -180.0:
                self.val_theta_4 = -180.0

            # revolute
            self.val_theta_5 = self.val_theta_5 + np.rad2deg(Joints_vel[4]) * 0.1
            if self.val_theta_5 > 180.0:
                self.val_theta_5 = 180.0
            elif self.val_theta_5 < -180.0:
                self.val_theta_5 = -180.0

            # revolute
            self.val_theta_6 = self.val_theta_6 + np.rad2deg(Joints_vel[5]) * 0.1
            if self.val_theta_6 > 180.0:
                self.val_theta_6 = 180.0
            elif self.val_theta_6 < -180.0:
                self.val_theta_6 = -180.0

            self.sld_theta_1.blockSignals(True)
            self.sld_theta_1.setValue(self.val_theta_1)
            self.sld_theta_1.blockSignals(False)

            self.sld_theta_2.blockSignals(True)
            self.sld_theta_2.setValue(self.val_theta_2)
            self.sld_theta_2.blockSignals(False)

            self.sld_theta_3.blockSignals(True)
            self.sld_theta_3.setValue(self.val_theta_3)
            self.sld_theta_3.blockSignals(False)

            self.sld_theta_4.blockSignals(True)
            self.sld_theta_4.setValue(self.val_theta_4)
            self.sld_theta_4.blockSignals(False)

            self.sld_theta_5.blockSignals(True)
            self.sld_theta_5.setValue(self.val_theta_5)
            self.sld_theta_5.blockSignals(False)

            self.sld_theta_6.blockSignals(True)
            self.sld_theta_6.setValue(self.val_theta_6)
            self.sld_theta_6.blockSignals(False)

        if self.scenario == 5:
            # Goal transformation (Body Frame-based)
            G_p_zyx = np.array([[self.val_pos_x], \
                                [self.val_pos_y], \
                                [self.val_pos_z]])
            G_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), G_p_zyx), axis=1), \
                                            np.array([[0, 0, 0, 1]])), axis=0)
            angle_z = np.deg2rad(self.val_rot_z)
            angle_y = np.deg2rad(self.val_rot_y)
            angle_x = np.deg2rad(self.val_rot_x)
            G_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0,                0,               1]])
            G_T_z = np.concatenate((np.concatenate((G_R_z, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            G_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            G_T_y = np.concatenate((np.concatenate((G_R_y, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            G_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
            G_T_x = np.concatenate((np.concatenate((G_R_x, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)

            G_T_zyx = (((self.S_T.dot(G_T_zyx_trans)).dot(G_T_z)).dot(G_T_y)).dot(G_T_x)

            G_T = G_T_zyx
            G_p = G_T_zyx[0:3,3]
            G_R = G_T_zyx[0:3,0:3]

            self.quiver_Gp = G_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Gx = G_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_Gy = G_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Gz = G_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_goal = 1.5
            self.axes.quiver(self.quiver_Gp[0], self.quiver_Gp[1], self.quiver_Gp[2], scale_goal*(self.quiver_Gx[0]-self.quiver_Gp[0]), scale_goal*(self.quiver_Gx[1]-self.quiver_Gp[1]), scale_goal*(self.quiver_Gx[2]-self.quiver_Gp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Gp[0], self.quiver_Gp[1], self.quiver_Gp[2], scale_goal*(self.quiver_Gy[0]-self.quiver_Gp[0]), scale_goal*(self.quiver_Gy[1]-self.quiver_Gp[1]), scale_goal*(self.quiver_Gy[2]-self.quiver_Gp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Gp[0], self.quiver_Gp[1], self.quiver_Gp[2], scale_goal*(self.quiver_Gz[0]-self.quiver_Gp[0]), scale_goal*(self.quiver_Gz[1]-self.quiver_Gp[1]), scale_goal*(self.quiver_Gz[2]-self.quiver_Gp[2]), color=['b'], arrow_length_ratio=0.15)

            # update goal for IK problem
            if self.fire_event:
                self.fire_event = False
 
                self.G_T_active = G_T

            if not np.isnan(self.G_T_active).any():
                # Solve IK problem at i-th iteration
                error_i = RigidBodyDynamics.Inverse(E_T).dot(self.G_T_active) 
                Twist_i = RigidBodyDynamics.LogMap(error_i)

                # Use Jacobian
                J_E_pinv = np.linalg.pinv(J_E) 
       
                Joints_vel = J_E_pinv.dot(Twist_i)

                # revolute
                self.val_theta_1 = self.val_theta_1 + np.rad2deg(Joints_vel[0]) * 0.1
                if self.val_theta_1 > 180.0:
                    self.val_theta_1 = 180.0
                elif self.val_theta_1 < -180.0:
                    self.val_theta_1 = -180.0

                # revolute
                self.val_theta_2 = self.val_theta_2 + np.rad2deg(Joints_vel[1]) * 0.1
                if self.val_theta_2 > 180.0:
                    self.val_theta_2 = 180.0
                elif self.val_theta_2 < -180.0:
                    self.val_theta_2 = -180.0

                # revolute
                self.val_theta_3 = self.val_theta_3 + np.rad2deg(Joints_vel[2]) * 0.1
                if self.val_theta_3 > 180.0:
                    self.val_theta_3 = 180.0
                elif self.val_theta_3 < -180.0:
                    self.val_theta_3 = -180.0

                # revolute
                self.val_theta_4 = self.val_theta_4 + np.rad2deg(Joints_vel[3]) * 0.1
                if self.val_theta_4 > 180.0:
                    self.val_theta_4 = 180.0
                elif self.val_theta_4 < -180.0:
                    self.val_theta_4 = -180.0

                # revolute
                self.val_theta_5 = self.val_theta_5 + np.rad2deg(Joints_vel[4]) * 0.1
                if self.val_theta_5 > 180.0:
                    self.val_theta_5 = 180.0
                elif self.val_theta_5 < -180.0:
                    self.val_theta_5 = -180.0

                # revolute
                self.val_theta_6 = self.val_theta_6 + np.rad2deg(Joints_vel[5]) * 0.1
                if self.val_theta_6 > 180.0:
                    self.val_theta_6 = 180.0
                elif self.val_theta_6 < -180.0:
                    self.val_theta_6 = -180.0

                self.sld_theta_1.blockSignals(True)
                self.sld_theta_1.setValue(self.val_theta_1)
                self.sld_theta_1.blockSignals(False)

                self.sld_theta_2.blockSignals(True)
                self.sld_theta_2.setValue(self.val_theta_2)
                self.sld_theta_2.blockSignals(False)

                self.sld_theta_3.blockSignals(True)
                self.sld_theta_3.setValue(self.val_theta_3)
                self.sld_theta_3.blockSignals(False)

                self.sld_theta_4.blockSignals(True)
                self.sld_theta_4.setValue(self.val_theta_4)
                self.sld_theta_4.blockSignals(False)

                self.sld_theta_5.blockSignals(True)
                self.sld_theta_5.setValue(self.val_theta_5)
                self.sld_theta_5.blockSignals(False)

                self.sld_theta_6.blockSignals(True)
                self.sld_theta_6.setValue(self.val_theta_6)
                self.sld_theta_6.blockSignals(False)

        self.canvas.draw()

    def on_update_values(self):       
        if self.cbx_theta_1.isChecked():
            self.cbx_theta_1.setChecked(False)
            self.sld_theta_1.setValue(0.0)
        if self.cbx_theta_2.isChecked():
            self.cbx_theta_2.setChecked(False)
            self.sld_theta_2.setValue(0.0)
        if self.cbx_theta_3.isChecked():
            self.cbx_theta_3.setChecked(False)
            self.sld_theta_3.setValue(0.0)
        if self.cbx_theta_4.isChecked():
            self.cbx_theta_4.setChecked(False)
            self.sld_theta_4.setValue(0.0)
        if self.cbx_theta_5.isChecked():
            self.cbx_theta_5.setChecked(False)
            self.sld_theta_5.setValue(0.0)
        if self.cbx_theta_6.isChecked():
            self.cbx_theta_6.setChecked(False)
            self.sld_theta_6.setValue(0.0)
        if self.cbx_vx.isChecked():
            self.cbx_vx.setChecked(False)
            self.sld_vx.setValue(0.0)
        if self.cbx_vy.isChecked():
            self.cbx_vy.setChecked(False)
            self.sld_vy.setValue(0.0)
        if self.cbx_vz.isChecked():
            self.cbx_vz.setChecked(False)
            self.sld_vz.setValue(0.0)
        if self.cbx_wx.isChecked():
            self.cbx_wx.setChecked(False)
            self.sld_wx.setValue(0.0)
        if self.cbx_wy.isChecked():
            self.cbx_wy.setChecked(False)
            self.sld_wy.setValue(0.0)
        if self.cbx_wz.isChecked():
            self.cbx_wz.setChecked(False)
            self.sld_wz.setValue(0.0)
        if self.cbx_posx.isChecked():
            self.cbx_posx.setChecked(False)
            self.sld_posx.setValue(0.0)
        if self.cbx_posy.isChecked():
            self.cbx_posy.setChecked(False)
            self.sld_posy.setValue(0.0)
        if self.cbx_posz.isChecked():
            self.cbx_posz.setChecked(False)
            self.sld_posz.setValue(0.0)
        if self.cbx_rotx.isChecked():
            self.cbx_rotx.setChecked(False)
            self.sld_rotx.setValue(0.0)
        if self.cbx_roty.isChecked():
            self.cbx_roty.setChecked(False)
            self.sld_roty.setValue(0.0)
        if self.cbx_rotz.isChecked():
            self.cbx_rotz.setChecked(False)
            self.sld_rotz.setValue(0.0)

        self.val_theta_1 = self.sld_theta_1.value()
        self.val_theta_2 = self.sld_theta_2.value()
        self.val_theta_3 = self.sld_theta_3.value()
        self.val_theta_4 = self.sld_theta_4.value()
        self.val_theta_5 = self.sld_theta_5.value()
        self.val_theta_6 = self.sld_theta_6.value()
        self.val_vel_x = self.sld_vx.value()
        self.val_vel_y = self.sld_vy.value()
        self.val_vel_z = self.sld_vz.value()
        self.val_omega_x = self.sld_wx.value()
        self.val_omega_y = self.sld_wy.value()
        self.val_omega_z = self.sld_wz.value()
        self.val_pos_x = self.sld_posx.value()
        self.val_pos_y = self.sld_posy.value()
        self.val_pos_z = self.sld_posz.value()
        self.val_rot_x = self.sld_rotx.value()
        self.val_rot_y = self.sld_roty.value()
        self.val_rot_z = self.sld_rotz.value()

        if self.rb_ik0.isChecked() and self.scenario != 0:
            self.scenario = 0
            self.rb_ik1.setChecked(False)
            self.rb_ik2.setChecked(False)
            self.rb_ik3.setChecked(False)
            self.rb_ik4.setChecked(False)
            self.rb_ik5.setChecked(False)

            self.fire_event = False
            self.G_T_active = np.NAN

            self.val_vel_x = 0
            self.val_vel_y = 0
            self.val_vel_z = 0
            self.val_omega_x = 0
            self.val_omega_y = 0
            self.val_omega_z = 0
            self.sld_vx.setValue(0)
            self.sld_vy.setValue(0)
            self.sld_vz.setValue(0)
            self.sld_wx.setValue(0)
            self.sld_wy.setValue(0)
            self.sld_wz.setValue(0)
            self.val_pos_x = 0
            self.val_pos_y = 0
            self.val_pos_z = 0
            self.val_rot_x = 0
            self.val_rot_y = 0
            self.val_rot_z = 0
            self.sld_posx.setValue(0)
            self.sld_posy.setValue(0)
            self.sld_posz.setValue(0)
            self.sld_rotx.setValue(0)
            self.sld_roty.setValue(0)
            self.sld_rotz.setValue(0)

            self.cbx_theta_1.setEnabled(True)
            self.sld_theta_1.setEnabled(True)
            self.cbx_theta_2.setEnabled(True)
            self.sld_theta_2.setEnabled(True)
            self.cbx_theta_3.setEnabled(True)
            self.sld_theta_3.setEnabled(True)
            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)
            self.cbx_vx.setEnabled(False)
            self.sld_vx.setEnabled(False)
            self.cbx_vy.setEnabled(False)
            self.sld_vy.setEnabled(False)
            self.cbx_vz.setEnabled(False)
            self.sld_vz.setEnabled(False)
            self.cbx_wx.setEnabled(False)
            self.sld_wx.setEnabled(False)
            self.cbx_wy.setEnabled(False)
            self.sld_wy.setEnabled(False)
            self.cbx_wz.setEnabled(False)
            self.sld_wz.setEnabled(False)
            self.cbx_posx.setEnabled(False)
            self.sld_posx.setEnabled(False)
            self.cbx_posy.setEnabled(False)
            self.sld_posy.setEnabled(False)
            self.cbx_posz.setEnabled(False)
            self.sld_posz.setEnabled(False)
            self.cbx_rotx.setEnabled(False)
            self.sld_rotx.setEnabled(False)
            self.cbx_roty.setEnabled(False)
            self.sld_roty.setEnabled(False)
            self.cbx_rotz.setEnabled(False)
            self.sld_rotz.setEnabled(False)
        elif self.rb_ik1.isChecked() and self.scenario != 1:
            self.scenario = 1
            self.rb_ik0.setChecked(False)
            self.rb_ik2.setChecked(False)
            self.rb_ik3.setChecked(False)
            self.rb_ik4.setChecked(False)
            self.rb_ik5.setChecked(False)

            self.fire_event = False
            self.G_T_active = np.NAN

            self.cbx_theta_1.setEnabled(False)
            self.sld_theta_1.setEnabled(False)
            self.cbx_theta_2.setEnabled(False)
            self.sld_theta_2.setEnabled(False)
            self.cbx_theta_3.setEnabled(False)
            self.sld_theta_3.setEnabled(False)
            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)
            self.cbx_vx.setEnabled(True)
            self.sld_vx.setEnabled(True)
            self.cbx_vy.setEnabled(True)
            self.sld_vy.setEnabled(True)
            self.cbx_vz.setEnabled(True)
            self.sld_vz.setEnabled(True)
            self.cbx_wx.setEnabled(True)
            self.sld_wx.setEnabled(True)
            self.cbx_wy.setEnabled(True)
            self.sld_wy.setEnabled(True)
            self.cbx_wz.setEnabled(True)
            self.sld_wz.setEnabled(True)
            self.cbx_posx.setEnabled(False)
            self.sld_posx.setEnabled(False)
            self.cbx_posy.setEnabled(False)
            self.sld_posy.setEnabled(False)
            self.cbx_posz.setEnabled(False)
            self.sld_posz.setEnabled(False)
            self.cbx_rotx.setEnabled(False)
            self.sld_rotx.setEnabled(False)
            self.cbx_roty.setEnabled(False)
            self.sld_roty.setEnabled(False)
            self.cbx_rotz.setEnabled(False)
            self.sld_rotz.setEnabled(False)
        elif self.rb_ik2.isChecked() and self.scenario != 2:
            self.scenario = 2
            self.rb_ik0.setChecked(False)
            self.rb_ik1.setChecked(False)
            self.rb_ik3.setChecked(False)
            self.rb_ik4.setChecked(False)
            self.rb_ik5.setChecked(False)

            self.fire_event = False
            self.G_T_active = np.NAN

            self.val_vel_x = 0
            self.val_vel_y = 0
            self.val_vel_z = 0
            self.val_omega_x = 0
            self.val_omega_y = 0
            self.val_omega_z = 0
            self.sld_vx.setValue(0)
            self.sld_vy.setValue(0)
            self.sld_vz.setValue(0)
            self.sld_wx.setValue(0)
            self.sld_wy.setValue(0)
            self.sld_wz.setValue(0)

            self.cbx_theta_1.setEnabled(False)
            self.sld_theta_1.setEnabled(False)
            self.cbx_theta_2.setEnabled(False)
            self.sld_theta_2.setEnabled(False)
            self.cbx_theta_3.setEnabled(False)
            self.sld_theta_3.setEnabled(False)
            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)
            self.cbx_vx.setEnabled(False)
            self.sld_vx.setEnabled(False)
            self.cbx_vy.setEnabled(False)
            self.sld_vy.setEnabled(False)
            self.cbx_vz.setEnabled(False)
            self.sld_vz.setEnabled(False)
            self.cbx_wx.setEnabled(False)
            self.sld_wx.setEnabled(False)
            self.cbx_wy.setEnabled(False)
            self.sld_wy.setEnabled(False)
            self.cbx_wz.setEnabled(False)
            self.sld_wz.setEnabled(False)
            self.cbx_posx.setEnabled(True)
            self.sld_posx.setEnabled(True)
            self.cbx_posy.setEnabled(True)
            self.sld_posy.setEnabled(True)
            self.cbx_posz.setEnabled(True)
            self.sld_posz.setEnabled(True)
            self.cbx_rotx.setEnabled(True)
            self.sld_rotx.setEnabled(True)
            self.cbx_roty.setEnabled(True)
            self.sld_roty.setEnabled(True)
            self.cbx_rotz.setEnabled(True)
            self.sld_rotz.setEnabled(True)

        elif self.rb_ik3.isChecked() and self.scenario != 3:
            self.scenario = 3
            self.rb_ik0.setChecked(False)
            self.rb_ik1.setChecked(False)
            self.rb_ik2.setChecked(False)
            self.rb_ik4.setChecked(False)
            self.rb_ik5.setChecked(False)

            self.fire_event = False
            self.G_T_active = np.NAN

            self.val_vel_x = 0
            self.val_vel_y = 0
            self.val_vel_z = 0
            self.val_omega_x = 0
            self.val_omega_y = 0
            self.val_omega_z = 0
            self.sld_vx.setValue(0)
            self.sld_vy.setValue(0)
            self.sld_vz.setValue(0)
            self.sld_wx.setValue(0)
            self.sld_wy.setValue(0)
            self.sld_wz.setValue(0)
            self.val_pos_x = 0
            self.val_pos_y = 0
            self.val_pos_z = 0
            self.val_rot_x = 0
            self.val_rot_y = 0
            self.val_rot_z = 0
            self.sld_posx.setValue(0)
            self.sld_posy.setValue(0)
            self.sld_posz.setValue(0)
            self.sld_rotx.setValue(0)
            self.sld_roty.setValue(0)
            self.sld_rotz.setValue(0)

            self.cbx_theta_1.setEnabled(True)
            self.sld_theta_1.setEnabled(True)
            self.cbx_theta_2.setEnabled(True)
            self.sld_theta_2.setEnabled(True)
            self.cbx_theta_3.setEnabled(True)
            self.sld_theta_3.setEnabled(True)
            self.cbx_theta_4.setEnabled(True)
            self.sld_theta_4.setEnabled(True)
            self.cbx_theta_5.setEnabled(True)
            self.sld_theta_5.setEnabled(True)
            self.cbx_theta_6.setEnabled(True)
            self.sld_theta_6.setEnabled(True)
            self.cbx_vx.setEnabled(False)
            self.sld_vx.setEnabled(False)
            self.cbx_vy.setEnabled(False)
            self.sld_vy.setEnabled(False)
            self.cbx_vz.setEnabled(False)
            self.sld_vz.setEnabled(False)
            self.cbx_wx.setEnabled(False)
            self.sld_wx.setEnabled(False)
            self.cbx_wy.setEnabled(False)
            self.sld_wy.setEnabled(False)
            self.cbx_wz.setEnabled(False)
            self.sld_wz.setEnabled(False)
            self.cbx_posx.setEnabled(False)
            self.sld_posx.setEnabled(False)
            self.cbx_posy.setEnabled(False)
            self.sld_posy.setEnabled(False)
            self.cbx_posz.setEnabled(False)
            self.sld_posz.setEnabled(False)
            self.cbx_rotx.setEnabled(False)
            self.sld_rotx.setEnabled(False)
            self.cbx_roty.setEnabled(False)
            self.sld_roty.setEnabled(False)
            self.cbx_rotz.setEnabled(False)
            self.sld_rotz.setEnabled(False)
        elif self.rb_ik4.isChecked() and self.scenario != 4:
            self.scenario = 4
            self.rb_ik0.setChecked(False)
            self.rb_ik1.setChecked(False)
            self.rb_ik2.setChecked(False)
            self.rb_ik3.setChecked(False)
            self.rb_ik5.setChecked(False)

            self.fire_event = False
            self.G_T_active = np.NAN

            self.cbx_theta_1.setEnabled(False)
            self.sld_theta_1.setEnabled(False)
            self.cbx_theta_2.setEnabled(False)
            self.sld_theta_2.setEnabled(False)
            self.cbx_theta_3.setEnabled(False)
            self.sld_theta_3.setEnabled(False)
            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)
            self.cbx_vx.setEnabled(True)
            self.sld_vx.setEnabled(True)
            self.cbx_vy.setEnabled(True)
            self.sld_vy.setEnabled(True)
            self.cbx_vz.setEnabled(True)
            self.sld_vz.setEnabled(True)
            self.cbx_wx.setEnabled(True)
            self.sld_wx.setEnabled(True)
            self.cbx_wy.setEnabled(True)
            self.sld_wy.setEnabled(True)
            self.cbx_wz.setEnabled(True)
            self.sld_wz.setEnabled(True)
            self.cbx_posx.setEnabled(False)
            self.sld_posx.setEnabled(False)
            self.cbx_posy.setEnabled(False)
            self.sld_posy.setEnabled(False)
            self.cbx_posz.setEnabled(False)
            self.sld_posz.setEnabled(False)
            self.cbx_rotx.setEnabled(False)
            self.sld_rotx.setEnabled(False)
            self.cbx_roty.setEnabled(False)
            self.sld_roty.setEnabled(False)
            self.cbx_rotz.setEnabled(False)
            self.sld_rotz.setEnabled(False)
        elif self.rb_ik5.isChecked() and self.scenario != 5:
            self.scenario = 5
            self.rb_ik0.setChecked(False)
            self.rb_ik1.setChecked(False)
            self.rb_ik2.setChecked(False)
            self.rb_ik3.setChecked(False)
            self.rb_ik4.setChecked(False)

            self.fire_event = False
            self.G_T_active = np.NAN

            self.val_vel_x = 0
            self.val_vel_y = 0
            self.val_vel_z = 0
            self.val_omega_x = 0
            self.val_omega_y = 0
            self.val_omega_z = 0
            self.sld_vx.setValue(0)
            self.sld_vy.setValue(0)
            self.sld_vz.setValue(0)
            self.sld_wx.setValue(0)
            self.sld_wy.setValue(0)
            self.sld_wz.setValue(0)

            self.cbx_theta_1.setEnabled(False)
            self.sld_theta_1.setEnabled(False)
            self.cbx_theta_2.setEnabled(False)
            self.sld_theta_2.setEnabled(False)
            self.cbx_theta_3.setEnabled(False)
            self.sld_theta_3.setEnabled(False)
            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)
            self.cbx_vx.setEnabled(False)
            self.sld_vx.setEnabled(False)
            self.cbx_vy.setEnabled(False)
            self.sld_vy.setEnabled(False)
            self.cbx_vz.setEnabled(False)
            self.sld_vz.setEnabled(False)
            self.cbx_wx.setEnabled(False)
            self.sld_wx.setEnabled(False)
            self.cbx_wy.setEnabled(False)
            self.sld_wy.setEnabled(False)
            self.cbx_wz.setEnabled(False)
            self.sld_wz.setEnabled(False)
            self.cbx_posx.setEnabled(True)
            self.sld_posx.setEnabled(True)
            self.cbx_posy.setEnabled(True)
            self.sld_posy.setEnabled(True)
            self.cbx_posz.setEnabled(True)
            self.sld_posz.setEnabled(True)
            self.cbx_rotx.setEnabled(True)
            self.sld_rotx.setEnabled(True)
            self.cbx_roty.setEnabled(True)
            self.sld_roty.setEnabled(True)
            self.cbx_rotz.setEnabled(True)
            self.sld_rotz.setEnabled(True)

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
        

        #self.draw_button = QPushButton("&Draw")
        #self.connect(self.draw_button, SIGNAL('clicked()'), self.on_update_values)
        
        self.cbx_theta_1 = QCheckBox('reset')
        self.cbx_theta_1.setChecked(False)
        self.connect(self.cbx_theta_1, SIGNAL('stateChanged(int)'), self.on_update_values)
        
        self.sld_theta_1 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_1.setMinimum(-180.0)
        self.sld_theta_1.setMaximum(180.0)
        self.sld_theta_1.setValue(0.0)
        self.sld_theta_1.setTracking(True)
        #self.sld_theta_1.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_1, SIGNAL('valueChanged(int)'), self.on_update_values)
        
        self.cbx_theta_2 = QCheckBox('reset')
        self.cbx_theta_2.setChecked(False)
        self.connect(self.cbx_theta_2, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_2 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_2.setMinimum(-180.0)
        self.sld_theta_2.setMaximum(180.0)
        self.sld_theta_2.setValue(0.0)
        self.sld_theta_2.setTracking(True)
        #self.sld_theta_2.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_2, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_3 = QCheckBox('reset')
        self.cbx_theta_3.setChecked(False)
        self.connect(self.cbx_theta_3, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_3 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_3.setMinimum(-180.0)
        self.sld_theta_3.setMaximum(180.0)
        self.sld_theta_3.setValue(0.0)
        self.sld_theta_3.setTracking(True)
        #self.sld_theta_3.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_3, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_4 = QCheckBox('reset')
        self.cbx_theta_4.setChecked(False)
        self.connect(self.cbx_theta_4, SIGNAL('stateChanged(int)'), self.on_update_values)
        
        self.sld_theta_4 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_4.setMinimum(-180.0)
        self.sld_theta_4.setMaximum(180.0)
        self.sld_theta_4.setValue(0.0)
        self.sld_theta_4.setTracking(True)
        #self.sld_theta_4.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_4, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_5 = QCheckBox('reset')
        self.cbx_theta_5.setChecked(False)
        self.connect(self.cbx_theta_5, SIGNAL('stateChanged(int)'), self.on_update_values)
        
        self.sld_theta_5 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_5.setMinimum(-180.0)
        self.sld_theta_5.setMaximum(180.0)
        self.sld_theta_5.setValue(0.0)
        self.sld_theta_5.setTracking(True)
        #self.sld_theta_5.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_5, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_6 = QCheckBox('reset')
        self.cbx_theta_6.setChecked(False)
        self.connect(self.cbx_theta_6, SIGNAL('stateChanged(int)'), self.on_update_values)
        
        self.sld_theta_6 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_6.setMinimum(-180.0)
        self.sld_theta_6.setMaximum(180.0)
        self.sld_theta_6.setValue(0.0)
        self.sld_theta_6.setTracking(True)
        #self.sld_theta_6.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_6, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_vx = QCheckBox('reset')
        self.cbx_vx.setChecked(False)
        self.connect(self.cbx_vx, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_vx = DoubleSlider(Qt.Horizontal)
        self.sld_vx.setMinimum(-1.0)
        self.sld_vx.setMaximum(1.0)
        self.sld_vx.setValue(0.0)
        self.sld_vx.setTracking(True)
        #self.sld_vx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vx, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_vy = QCheckBox('reset')
        self.cbx_vy.setChecked(False)
        self.connect(self.cbx_vy, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_vy = DoubleSlider(Qt.Horizontal)
        self.sld_vy.setMinimum(-1.0)
        self.sld_vy.setMaximum(1.0)
        self.sld_vy.setValue(0.0)
        self.sld_vy.setTracking(True)
        #self.sld_vy.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vy, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_vz = QCheckBox('reset')
        self.cbx_vz.setChecked(False)
        self.connect(self.cbx_vz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_vz = DoubleSlider(Qt.Horizontal)
        self.sld_vz.setMinimum(-1.0)
        self.sld_vz.setMaximum(1.0)
        self.sld_vz.setValue(0.0)
        self.sld_vz.setTracking(True)
        #self.sld_vz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_wx = QCheckBox('reset')
        self.cbx_wx.setChecked(False)
        self.connect(self.cbx_wx, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_wx = DoubleSlider(Qt.Horizontal)
        self.sld_wx.setMinimum(-180.0)
        self.sld_wx.setMaximum(180.0)
        self.sld_wx.setValue(0.0)
        self.sld_wx.setTracking(True)
        #self.sld_wx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_wx, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_wy = QCheckBox('reset')
        self.cbx_wy.setChecked(False)
        self.connect(self.cbx_wy, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_wy = DoubleSlider(Qt.Horizontal)
        self.sld_wy.setMinimum(-180.0)
        self.sld_wy.setMaximum(180.0)
        self.sld_wy.setValue(0.0)
        self.sld_wy.setTracking(True)
        #self.sld_wy.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_wy, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_wz = QCheckBox('reset')
        self.cbx_wz.setChecked(False)
        self.connect(self.cbx_wz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_wz = DoubleSlider(Qt.Horizontal)
        self.sld_wz.setMinimum(-180.0)
        self.sld_wz.setMaximum(180.0)
        self.sld_wz.setValue(0.0)
        self.sld_wz.setTracking(True)
        #self.sld_wz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_wz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_posx = QCheckBox('reset')
        self.cbx_posx.setChecked(False)
        self.connect(self.cbx_posx, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_posx = DoubleSlider(Qt.Horizontal)
        self.sld_posx.setMinimum(-2.0)
        self.sld_posx.setMaximum(2.0)
        self.sld_posx.setValue(0.0)
        self.sld_posx.setTracking(True)
        #self.sld_posx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_posx, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_posy = QCheckBox('reset')
        self.cbx_posy.setChecked(False)
        self.connect(self.cbx_posy, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_posy = DoubleSlider(Qt.Horizontal)
        self.sld_posy.setMinimum(-2.0)
        self.sld_posy.setMaximum(2.0)
        self.sld_posy.setValue(0.0)
        self.sld_posy.setTracking(True)
        #self.sld_posy.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_posy, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_posz = QCheckBox('reset')
        self.cbx_posz.setChecked(False)
        self.connect(self.cbx_posz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_posz = DoubleSlider(Qt.Horizontal)
        self.sld_posz.setMinimum(-2.0)
        self.sld_posz.setMaximum(2.0)
        self.sld_posz.setValue(0.0)
        self.sld_posz.setTracking(True)
        #self.sld_posz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_posz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_rotx = QCheckBox('reset')
        self.cbx_rotx.setChecked(False)
        self.connect(self.cbx_rotx, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_rotx = DoubleSlider(Qt.Horizontal)
        self.sld_rotx.setMinimum(-180.0)
        self.sld_rotx.setMaximum(180.0)
        self.sld_rotx.setValue(0.0)
        self.sld_rotx.setTracking(True)
        #self.sld_rotx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_rotx, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_roty = QCheckBox('reset')
        self.cbx_roty.setChecked(False)
        self.connect(self.cbx_roty, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_roty = DoubleSlider(Qt.Horizontal)
        self.sld_roty.setMinimum(-180.0)
        self.sld_roty.setMaximum(180.0)
        self.sld_roty.setValue(0.0)
        self.sld_roty.setTracking(True)
        #self.sld_roty.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_roty, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_rotz = QCheckBox('reset')
        self.cbx_rotz.setChecked(False)
        self.connect(self.cbx_rotz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_rotz = DoubleSlider(Qt.Horizontal)
        self.sld_rotz.setMinimum(-180.0)
        self.sld_rotz.setMaximum(180.0)
        self.sld_rotz.setValue(0.0)
        self.sld_rotz.setTracking(True)
        #self.sld_rotz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_rotz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.rb_ik0 = QCheckBox('FK #1')
        self.rb_ik0.setChecked(True)
        self.connect(self.rb_ik0, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_ik1 = QCheckBox('IK #1')
        self.rb_ik1.setChecked(False)
        self.connect(self.rb_ik1, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_ik2 = QCheckBox('Full IK #1')
        self.rb_ik2.setChecked(False)
        self.connect(self.rb_ik2, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_ik3 = QCheckBox('FK #2')
        self.rb_ik3.setChecked(False)
        self.connect(self.rb_ik3, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_ik4 = QCheckBox('IK #2')
        self.rb_ik4.setChecked(False)
        self.connect(self.rb_ik4, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_ik5 = QCheckBox('Full IK #2')
        self.rb_ik5.setChecked(False)
        self.connect(self.rb_ik5, SIGNAL('stateChanged(int)'), self.on_update_values)

        hbox_theta1 = QHBoxLayout()
        for w in [ self.cbx_theta_1, QLabel('θ1'), QLabel('-180'), self.sld_theta_1, QLabel('180')]:
            hbox_theta1.addWidget(w)
            hbox_theta1.setAlignment(w, Qt.AlignVCenter)

        hbox_theta2 = QHBoxLayout()
        for w in [ self.cbx_theta_2, QLabel('θ2'), QLabel('-180'), self.sld_theta_2, QLabel('180')]:
            hbox_theta2.addWidget(w)
            hbox_theta2.setAlignment(w, Qt.AlignVCenter)

        hbox_theta3 = QHBoxLayout()
        for w in [ self.cbx_theta_3, QLabel('θ3'), QLabel('-180'), self.sld_theta_3, QLabel('180')]:
            hbox_theta3.addWidget(w)
            hbox_theta3.setAlignment(w, Qt.AlignVCenter)

        hbox_theta4 = QHBoxLayout()
        for w in [ self.cbx_theta_4, QLabel('θ4'), QLabel('-180'), self.sld_theta_4, QLabel('180')]:
            hbox_theta4.addWidget(w)
            hbox_theta4.setAlignment(w, Qt.AlignVCenter)

        hbox_theta5 = QHBoxLayout()
        for w in [ self.cbx_theta_5, QLabel('θ5'), QLabel('-180'), self.sld_theta_5, QLabel('180')]:
            hbox_theta5.addWidget(w)
            hbox_theta5.setAlignment(w, Qt.AlignVCenter)

        hbox_theta6 = QHBoxLayout()
        for w in [ self.cbx_theta_6, QLabel('θ6'), QLabel('-180'), self.sld_theta_6, QLabel('180')]:
            hbox_theta6.addWidget(w)
            hbox_theta6.setAlignment(w, Qt.AlignVCenter)

        hbox_vx = QHBoxLayout()
        for w in [ self.cbx_vx, QLabel('vx'), QLabel('-1'), self.sld_vx, QLabel('1'), self.cbx_posx, QLabel('px'), QLabel('-2'), self.sld_posx, QLabel('2')]:
            hbox_vx.addWidget(w)
            hbox_vx.setAlignment(w, Qt.AlignVCenter)

        hbox_vy = QHBoxLayout()
        for w in [ self.cbx_vy, QLabel('vy'), QLabel('-1'), self.sld_vy, QLabel('1'), self.cbx_posy, QLabel('py'), QLabel('-2'), self.sld_posy, QLabel('2')]:
            hbox_vy.addWidget(w)
            hbox_vy.setAlignment(w, Qt.AlignVCenter)

        hbox_vz = QHBoxLayout()
        for w in [ self.cbx_vz, QLabel('vz'), QLabel('-1'), self.sld_vz, QLabel('1'), self.cbx_posz, QLabel('pz'), QLabel('-2'), self.sld_posz, QLabel('2')]:
            hbox_vz.addWidget(w)
            hbox_vz.setAlignment(w, Qt.AlignVCenter)

        hbox_wx = QHBoxLayout()
        for w in [ self.cbx_wx, QLabel('ωx'), QLabel('-180'), self.sld_wx, QLabel('180'), self.cbx_rotx, QLabel('rx'), QLabel('-180'), self.sld_rotx, QLabel('180')]:
            hbox_wx.addWidget(w)
            hbox_wx.setAlignment(w, Qt.AlignVCenter)

        hbox_wy = QHBoxLayout()
        for w in [ self.cbx_wy, QLabel('ωy'), QLabel('-180'), self.sld_wy, QLabel('180'), self.cbx_roty, QLabel('ry'), QLabel('-180'), self.sld_roty, QLabel('180')]:
            hbox_wy.addWidget(w)
            hbox_wy.setAlignment(w, Qt.AlignVCenter)

        hbox_wz = QHBoxLayout()
        for w in [ self.cbx_wz, QLabel('ωz'), QLabel('-180'), self.sld_wz, QLabel('180'), self.cbx_rotz, QLabel('rz'), QLabel('-180'), self.sld_rotz, QLabel('180')]:
            hbox_wz.addWidget(w)
            hbox_wz.setAlignment(w, Qt.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ self.rb_ik0, self.rb_ik1, self.rb_ik2, self.rb_ik3, self.rb_ik4, self.rb_ik5 ]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_theta1)
        vbox.addLayout(hbox_theta2)
        vbox.addLayout(hbox_theta3)
        vbox.addLayout(hbox_theta4)
        vbox.addLayout(hbox_theta5)
        vbox.addLayout(hbox_theta6)
        vbox.addLayout(hbox_vx)
        vbox.addLayout(hbox_vy)
        vbox.addLayout(hbox_vz)
        vbox.addLayout(hbox_wx)
        vbox.addLayout(hbox_wy)
        vbox.addLayout(hbox_wz)
        vbox.addLayout(hbox_rb)

        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Inverse Kinematics")
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

