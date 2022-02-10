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
        self.setWindowTitle('Forward Kinematics')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()


        self.val_theta_1 = 0
        self.val_theta_2 = 0
        self.val_theta_3 = 0
        self.val_theta_4 = 0
        self.val_theta_5 = 0
        self.val_theta_6 = 0

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

        if self.scenario == 0:
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

            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0       , -S1[2][0],  S1[1][0]], \
                                      [ S1[2][0],         0, -S1[0][0]], \
                                      [-S1[1][0],  S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_2)/np.pi) < 1e-3:
                e_S2Matrix_theta2 = np.eye(4)
            elif np.linalg.norm(S2[0:3]) < 1e-6:
                V = np.eye(3);
                e_S2Matrix_theta2 = np.concatenate((np.concatenate((np.eye(3), V.dot(S2[3:6]) * np.deg2rad(self.val_theta_2)/np.pi), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0       , -S2[2][0],  S2[1][0]], \
                                      [ S2[2][0],         0, -S2[0][0]], \
                                      [-S2[1][0],  S2[0][0],         0]]) * np.deg2rad(self.val_theta_2)/np.pi
                theta = np.linalg.norm(S2[0:3] * np.deg2rad(self.val_theta_2)/np.pi)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_S2Matrix_theta2 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(S2[3:6] * np.deg2rad(self.val_theta_2)/np.pi)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_3)) < 1e-3:
                e_S3Matrix_theta3 = np.eye(4)
            elif np.linalg.norm(S3[0:3]) < 1e-6:
                V = np.eye(3);
                e_S3Matrix_theta3 = np.concatenate((np.concatenate((np.eye(3), V.dot(S3[3:6]) * np.deg2rad(self.val_theta_3)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0       , -S3[2][0],  S3[1][0]], \
                                      [ S3[2][0],         0, -S3[0][0]], \
                                      [-S3[1][0],  S3[0][0],         0]]) * np.deg2rad(self.val_theta_3)
                theta = np.linalg.norm(S3[0:3] * np.deg2rad(self.val_theta_3))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_S3Matrix_theta3 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(S3[3:6] * np.deg2rad(self.val_theta_3))), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)

            E_T = ((Home_T.dot(e_S1Matrix_theta1)).dot(e_S2Matrix_theta2)).dot(e_S3Matrix_theta3)  # Post-multiply

            # Intermediate frames
            # Joint 1
            Home_j1 = np.concatenate((np.concatenate((np.eye(3), np.array([[(0)], [0], [0]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)
             
            j1_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(0)], [0], [0]]), axis=0)), axis=0)
           
            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_j1S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(j1_S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_j1S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(j1_S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j1_S1[2][0],  j1_S1[1][0]], \
                                      [ j1_S1[2][0],            0, -j1_S1[0][0]], \
                                      [-j1_S1[1][0],  j1_S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(j1_S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j1S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j1_S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

            j1_T = Home_j1.dot(e_j1S1Matrix_theta1)  # Post-multiply

            # Joint 2
            Home_j2 = np.concatenate((np.concatenate((np.eye(3), np.array([[(L1+L2)], [0], [0]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j2_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(L1+L2)], [0], [0]]), axis=0)), axis=0)
            j2_S2 = np.concatenate((np.array([[0], [0], [0]]), \
                                    np.array([[1], [0], [0]])), axis=0)
           
            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_j2S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(j2_S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_j2S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S1[2][0],  j2_S1[1][0]], \
                                      [ j2_S1[2][0],            0, -j2_S1[0][0]], \
                                      [-j2_S1[1][0],  j2_S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(j2_S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j2S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_2)/np.pi) < 1e-3:
                e_j2S2Matrix_theta2 = np.eye(4)
            elif np.linalg.norm(j2_S2[0:3]) < 1e-6:
                V = np.eye(3);
                e_j2S2Matrix_theta2 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S2[3:6]) * np.deg2rad(self.val_theta_2)/np.pi), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S2[2][0],  j2_S2[1][0]], \
                                      [ j2_S2[2][0],            0, -j2_S2[0][0]], \
                                      [-j2_S2[1][0],  j2_S2[0][0],         0]]) * np.deg2rad(self.val_theta_2)/np.pi
                theta = np.linalg.norm(j2_S2[0:3] * np.deg2rad(self.val_theta_2)/np.pi)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j2S2Matrix_theta2 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S2[3:6] * np.deg2rad(self.val_theta_2)/np.pi)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
  
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

            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_j3S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(j2_S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_j3S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S1[2][0],  j2_S1[1][0]], \
                                      [ j2_S1[2][0],            0, -j2_S1[0][0]], \
                                      [-j2_S1[1][0],  j2_S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(j2_S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j3S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_2)/np.pi) < 1e-3:
                e_j3S2Matrix_theta2 = np.eye(4)
            elif np.linalg.norm(j2_S2[0:3]) < 1e-6:
                V = np.eye(3);
                e_j3S2Matrix_theta2 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S2[3:6]) * np.deg2rad(self.val_theta_2)/np.pi), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S2[2][0],  j2_S2[1][0]], \
                                      [ j2_S2[2][0],            0, -j2_S2[0][0]], \
                                      [-j2_S2[1][0],  j2_S2[0][0],         0]]) * np.deg2rad(self.val_theta_2)/np.pi
                theta = np.linalg.norm(j2_S2[0:3] * np.deg2rad(self.val_theta_2)/np.pi)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j3S2Matrix_theta2 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S2[3:6] * np.deg2rad(self.val_theta_2)/np.pi)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
  
            if abs(np.deg2rad(self.val_theta_3)) < 1e-3:
                e_j3S3Matrix_theta3 = np.eye(4)
            elif np.linalg.norm(j3_S3[0:3]) < 1e-6:
                V = np.eye(3);
                e_j3S3Matrix_theta3 = np.concatenate((np.concatenate((np.eye(3), V.dot(j3_S3[3:6]) * np.deg2rad(self.val_theta_3)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j3_S3[2][0],  j3_S3[1][0]], \
                                      [ j3_S3[2][0],            0, -j3_S3[0][0]], \
                                      [-j3_S3[1][0],  j3_S3[0][0],         0]]) * np.deg2rad(self.val_theta_3)
                theta = np.linalg.norm(j3_S3[0:3] * np.deg2rad(self.val_theta_3))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j3S3Matrix_theta3 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j3_S3[3:6] * np.deg2rad(self.val_theta_3))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

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

        if self.scenario == 1:
            # 3-Joint SE(3)
            L1 = 0.25
            L2 = 0.10
            L3 = 0.75
            L4 = 0.10
            L5 = 0.50

            Home_T = np.concatenate((np.concatenate((np.eye(3), np.array([[(-L4)], [(L2)], [(L1+L3+L5)]])), axis=1), \
                                     np.array([[0, 0, 0, 1]])), axis=0)

            S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                 np.cross(-np.array([[0], [0], [1]]) , np.array([[L4], [-(L2)], [-(L1+L3+L5)]]), axis=0)), axis=0)
            S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                 np.cross(-np.array([[0], [1], [0]]) , np.array([[L4], [-(0)], [-(L3+L5)]]), axis=0)), axis=0)
            S3 = np.concatenate((np.array([[1], [0], [0]]), \
                                 np.cross(-np.array([[1], [0], [0]]) , np.array([[0], [-(0)], [-(L5)]]), axis=0)), axis=0)

            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0       , -S1[2][0],  S1[1][0]], \
                                      [ S1[2][0],         0, -S1[0][0]], \
                                      [-S1[1][0],  S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_2)/np.pi) < 1e-3:
                e_S2Matrix_theta2 = np.eye(4)
            elif np.linalg.norm(S2[0:3]) < 1e-6:
                V = np.eye(3);
                e_S2Matrix_theta2 = np.concatenate((np.concatenate((np.eye(3), V.dot(S2[3:6]) * np.deg2rad(self.val_theta_2)/np.pi), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0       , -S2[2][0],  S2[1][0]], \
                                      [ S2[2][0],         0, -S2[0][0]], \
                                      [-S2[1][0],  S2[0][0],         0]]) * np.deg2rad(self.val_theta_2)/np.pi
                theta = np.linalg.norm(S2[0:3] * np.deg2rad(self.val_theta_2)/np.pi)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_S2Matrix_theta2 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(S2[3:6] * np.deg2rad(self.val_theta_2)/np.pi)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_3)) < 1e-3:
                e_S3Matrix_theta3 = np.eye(4)
            elif np.linalg.norm(S3[0:3]) < 1e-6:
                V = np.eye(3);
                e_S3Matrix_theta3 = np.concatenate((np.concatenate((np.eye(3), V.dot(S3[3:6]) * np.deg2rad(self.val_theta_3)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0       , -S3[2][0],  S3[1][0]], \
                                      [ S3[2][0],         0, -S3[0][0]], \
                                      [-S3[1][0],  S3[0][0],         0]]) * np.deg2rad(self.val_theta_3)
                theta = np.linalg.norm(S3[0:3] * np.deg2rad(self.val_theta_3))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_S3Matrix_theta3 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(S3[3:6] * np.deg2rad(self.val_theta_3))), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)

            E_T = ((Home_T.dot(e_S1Matrix_theta1)).dot(e_S2Matrix_theta2)).dot(e_S3Matrix_theta3)  # Post-multiply

            # Intermediate frames
            # Joint 1
            Home_j1 = np.concatenate((np.concatenate((np.eye(3), np.array([[-(0)], [(0)], [(0)]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)
             
            j1_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(0)], [0], [0]]), axis=0)), axis=0)
           
            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_j1S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(j1_S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_j1S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(j1_S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j1_S1[2][0],  j1_S1[1][0]], \
                                      [ j1_S1[2][0],            0, -j1_S1[0][0]], \
                                      [-j1_S1[1][0],  j1_S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(j1_S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j1S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j1_S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

            j1_T = Home_j1.dot(e_j1S1Matrix_theta1)  # Post-multiply

            # Joint 2
            Home_j2 = np.concatenate((np.concatenate((np.eye(3), np.array([[-(0)], [(L2)], [(L1)]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j2_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[-(0)], [-(L2)], [-(L1)]]), axis=0)), axis=0)
            j2_S2 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[0], [-(0)], [-(0)]]), axis=0)), axis=0)
           
            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_j2S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(j2_S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_j2S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S1[2][0],  j2_S1[1][0]], \
                                      [ j2_S1[2][0],            0, -j2_S1[0][0]], \
                                      [-j2_S1[1][0],  j2_S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(j2_S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j2S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_2)/np.pi) < 1e-3:
                e_j2S2Matrix_theta2 = np.eye(4)
            elif np.linalg.norm(j2_S2[0:3]) < 1e-6:
                V = np.eye(3);
                e_j2S2Matrix_theta2 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S2[3:6]) * np.deg2rad(self.val_theta_2)/np.pi), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S2[2][0],  j2_S2[1][0]], \
                                      [ j2_S2[2][0],            0, -j2_S2[0][0]], \
                                      [-j2_S2[1][0],  j2_S2[0][0],         0]]) * np.deg2rad(self.val_theta_2)/np.pi
                theta = np.linalg.norm(j2_S2[0:3] * np.deg2rad(self.val_theta_2)/np.pi)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j2S2Matrix_theta2 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S2[3:6] * np.deg2rad(self.val_theta_2)/np.pi)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
  
            j2_T = (Home_j2.dot(e_j2S1Matrix_theta1)).dot(e_j2S2Matrix_theta2)  # Post-multiply
            
            # Joint 3
            Home_j3 = np.concatenate((np.concatenate((np.eye(3), np.array([[-(L4)], [(L2)], [(L1+L3)]])), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            j3_S1 = np.concatenate((np.array([[0], [0], [1]]), \
                                    np.cross(-np.array([[0], [0], [1]]) , np.array([[L4], [-(L2)], [-(L1+L3)]]), axis=0)), axis=0)
            j3_S1 = np.concatenate((np.array([[0], [1], [0]]), \
                                    np.cross(-np.array([[0], [1], [0]]) , np.array([[L4], [-(0)], [-(L3)]]), axis=0)), axis=0)
            j3_S3 = np.concatenate((np.array([[1], [0], [0]]), \
                                    np.cross(-np.array([[1], [0], [0]]) , np.array([[0], [-(0)], [-(0)]]), axis=0)), axis=0)     

            if abs(np.deg2rad(self.val_theta_1)) < 1e-3:
                e_j3S1Matrix_theta1 = np.eye(4)
            elif np.linalg.norm(j2_S1[0:3]) < 1e-6:
                V = np.eye(3);
                e_j3S1Matrix_theta1 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S1[3:6]) * np.deg2rad(self.val_theta_1)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S1[2][0],  j2_S1[1][0]], \
                                      [ j2_S1[2][0],            0, -j2_S1[0][0]], \
                                      [-j2_S1[1][0],  j2_S1[0][0],         0]]) * np.deg2rad(self.val_theta_1)
                theta = np.linalg.norm(j2_S1[0:3] * np.deg2rad(self.val_theta_1))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j3S1Matrix_theta1 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S1[3:6] * np.deg2rad(self.val_theta_1))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

            if abs(np.deg2rad(self.val_theta_2)/np.pi) < 1e-3:
                e_j3S2Matrix_theta2 = np.eye(4)
            elif np.linalg.norm(j2_S2[0:3]) < 1e-6:
                V = np.eye(3);
                e_j3S2Matrix_theta2 = np.concatenate((np.concatenate((np.eye(3), V.dot(j2_S2[3:6]) * np.deg2rad(self.val_theta_2)/np.pi), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j2_S2[2][0],  j2_S2[1][0]], \
                                      [ j2_S2[2][0],            0, -j2_S2[0][0]], \
                                      [-j2_S2[1][0],  j2_S2[0][0],         0]]) * np.deg2rad(self.val_theta_2)/np.pi
                theta = np.linalg.norm(j2_S2[0:3] * np.deg2rad(self.val_theta_2)/np.pi)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j3S2Matrix_theta2 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j2_S2[3:6] * np.deg2rad(self.val_theta_2)/np.pi)), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)
  
            if abs(np.deg2rad(self.val_theta_3)) < 1e-3:
                e_j3S3Matrix_theta3 = np.eye(4)
            elif np.linalg.norm(j3_S3[0:3]) < 1e-6:
                V = np.eye(3);
                e_j3S3Matrix_theta3 = np.concatenate((np.concatenate((np.eye(3), V.dot(j3_S3[3:6]) * np.deg2rad(self.val_theta_3)), axis=1), \
                                                  np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[ 0          , -j3_S3[2][0],  j3_S3[1][0]], \
                                      [ j3_S3[2][0],            0, -j3_S3[0][0]], \
                                      [-j3_S3[1][0],  j3_S3[0][0],         0]]) * np.deg2rad(self.val_theta_3)
                theta = np.linalg.norm(j3_S3[0:3] * np.deg2rad(self.val_theta_3))
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_j3S3Matrix_theta3 = np.concatenate((np.concatenate((e_omegaSkew, V.dot(j3_S3[3:6] * np.deg2rad(self.val_theta_3))), axis=1), \
                                                      np.array([[0, 0, 0, 1]])), axis=0)

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

        self.val_theta_1 = self.sld_theta_1.value()
        self.val_theta_2 = self.sld_theta_2.value()
        self.val_theta_3 = self.sld_theta_3.value()
        self.val_theta_4 = self.sld_theta_4.value()
        self.val_theta_5 = self.sld_theta_5.value()
        self.val_theta_6 = self.sld_theta_6.value()

        if self.rb_fk0.isChecked() and self.scenario != 0:
            self.scenario = 0
            self.rb_fk1.setChecked(False)
            
            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)
        elif self.rb_fk1.isChecked() and self.scenario != 1:
            self.scenario = 1
            self.rb_fk0.setChecked(False)

            self.cbx_theta_4.setEnabled(False)
            self.sld_theta_4.setEnabled(False)
            self.cbx_theta_5.setEnabled(False)
            self.sld_theta_5.setEnabled(False)
            self.cbx_theta_6.setEnabled(False)
            self.sld_theta_6.setEnabled(False)

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
        self.sld_theta_1.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_1, SIGNAL('valueChanged(int)'), self.on_update_values)
        
        self.cbx_theta_2 = QCheckBox('reset')
        self.cbx_theta_2.setChecked(False)
        self.connect(self.cbx_theta_2, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_2 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_2.setMinimum(-180.0)
        self.sld_theta_2.setMaximum(180.0)
        self.sld_theta_2.setValue(0.0)
        self.sld_theta_2.setTracking(True)
        self.sld_theta_2.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_2, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_3 = QCheckBox('reset')
        self.cbx_theta_3.setChecked(False)
        self.connect(self.cbx_theta_3, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_3 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_3.setMinimum(-180.0)
        self.sld_theta_3.setMaximum(180.0)
        self.sld_theta_3.setValue(0.0)
        self.sld_theta_3.setTracking(True)
        self.sld_theta_3.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_3, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_4 = QCheckBox('reset')
        self.cbx_theta_4.setChecked(False)
        self.connect(self.cbx_theta_4, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_4 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_4.setMinimum(-1.0)
        self.sld_theta_4.setMaximum(1.0)
        self.sld_theta_4.setValue(0.0)
        self.sld_theta_4.setTracking(True)
        self.sld_theta_4.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_4, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_5 = QCheckBox('reset')
        self.cbx_theta_5.setChecked(False)
        self.connect(self.cbx_theta_5, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_5 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_5.setMinimum(-1.0)
        self.sld_theta_5.setMaximum(1.0)
        self.sld_theta_5.setValue(0.0)
        self.sld_theta_5.setTracking(True)
        self.sld_theta_5.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_5, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_theta_6 = QCheckBox('reset')
        self.cbx_theta_6.setChecked(False)
        self.connect(self.cbx_theta_6, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_theta_6 = DoubleSlider(Qt.Horizontal)
        self.sld_theta_6.setMinimum(-1.0)
        self.sld_theta_6.setMaximum(1.0)
        self.sld_theta_6.setValue(0.0)
        self.sld_theta_6.setTracking(True)
        self.sld_theta_6.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_theta_6, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.rb_fk0 = QCheckBox('FK #1')
        self.rb_fk0.setChecked(True)
        self.connect(self.rb_fk0, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_fk1 = QCheckBox('FK #2')
        self.rb_fk1.setChecked(False)
        self.connect(self.rb_fk1, SIGNAL('stateChanged(int)'), self.on_update_values)

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
        for w in [ self.cbx_theta_4, QLabel('θ4'), QLabel('-1'), self.sld_theta_4, QLabel('1')]:
            hbox_theta4.addWidget(w)
            hbox_theta4.setAlignment(w, Qt.AlignVCenter)

        hbox_theta5 = QHBoxLayout()
        for w in [ self.cbx_theta_5, QLabel('θ5'), QLabel('-1'), self.sld_theta_5, QLabel('1')]:
            hbox_theta5.addWidget(w)
            hbox_theta5.setAlignment(w, Qt.AlignVCenter)

        hbox_theta6 = QHBoxLayout()
        for w in [ self.cbx_theta_6, QLabel('θ6'), QLabel('-1'), self.sld_theta_6, QLabel('1')]:
            hbox_theta6.addWidget(w)
            hbox_theta6.setAlignment(w, Qt.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ self.rb_fk0, self.rb_fk1 ]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_theta1)
        vbox.addLayout(hbox_theta2)
        vbox.addLayout(hbox_theta3)
        vbox.addLayout(hbox_theta4)
        vbox.addLayout(hbox_theta5)
        vbox.addLayout(hbox_theta6)
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

