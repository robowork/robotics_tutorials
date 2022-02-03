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
        self.setWindowTitle('Rigid Body Motion: Twists & Kinematics')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()


        self.val_omega_x = 0
        self.val_omega_y = 0
        self.val_omega_z = 0
        self.val_vel_x = 0
        self.val_vel_y = 0
        self.val_vel_z = 0

        self.B_p_x = 0.50
        self.B_p_y = 0.25
        self.B_p_z = 0.75
        self.B_p_qw = 0.872
        self.B_p_qx = 0.215
        self.B_p_qy = 0.189
        self.B_p_qz = 0.398

        self.ef_tx.setText(str(self.B_p_x))
        self.ef_ty.setText(str(self.B_p_y))
        self.ef_tz.setText(str(self.B_p_z))
        self.ef_qw.setText(str(self.B_p_qw))
        self.ef_qx.setText(str(self.B_p_qx))
        self.ef_qy.setText(str(self.B_p_qy))
        self.ef_qz.setText(str(self.B_p_qz))

        self.timer_period = 0.10
        self.timer_value = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_timer)
        self.timer.start(self.timer_period * 1000)

        self.scenario = 0
        
        # Space Frame / origin
        self.S_p = np.array([[0], [0], [0]])
        self.S_R = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.S_T = np.concatenate((np.concatenate((self.S_R, self.S_p), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
   
        self.quiver_Sx = np.array([[1],[0],[0]])
        self.quiver_Sy = np.array([[0],[1],[0]])
        self.quiver_Sz = np.array([[0],[0],[1]])
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
            # Twist expressed w.r.t. Body Frame
            B_p_zyx = np.array([[self.B_p_x], [self.B_p_y], [self.B_p_z]])
            B_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), B_p_zyx), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx_rot = np.concatenate((np.concatenate((quaternions.quat2mat(np.array([self.B_p_qw, self.B_p_qx, self.B_p_qy, self.B_p_qz])), np.zeros((3,1))), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx = (self.S_T.dot(B_T_zyx_trans)).dot(B_T_zyx_rot)

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            Twist_B = np.array([[np.deg2rad(self.val_omega_x)], [np.deg2rad(self.val_omega_y)], [np.deg2rad(self.val_omega_z)], [self.val_vel_x], [self.val_vel_y], [self.val_vel_z]])

            if np.linalg.norm(Twist_B[0:3]) < 1e-6:
                V = np.eye(3)
                e_VMatrix_theta = np.concatenate((np.concatenate((np.eye(3), V.dot(Twist_B[3:6]) * self.timer_value), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[0, -Twist_B[2][0], Twist_B[1][0]], [Twist_B[2][0], 0, -Twist_B[0][0]], [-Twist_B[1][0], Twist_B[0][0], 0]]) * self.timer_value
                theta = np.linalg.norm(Twist_B[0:3] * self.timer_value)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_VMatrix_theta = np.concatenate((np.concatenate((e_omegaSkew, V.dot(Twist_B[3:6] * self.timer_value)), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
      
            B_T = B_T_zyx.dot(e_VMatrix_theta)  # Post-multiply

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

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

        elif self.scenario == 1:
            # Twist expressed w.r.t. Space Frame
            B_p_zyx = np.array([[self.B_p_x], [self.B_p_y], [self.B_p_z]])
            B_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), B_p_zyx), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx_rot = np.concatenate((np.concatenate((quaternions.quat2mat(np.array([self.B_p_qw, self.B_p_qx, self.B_p_qy, self.B_p_qz])), np.zeros((3,1))), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx = (self.S_T.dot(B_T_zyx_trans)).dot(B_T_zyx_rot)

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            Twist_B = np.array([[np.deg2rad(self.val_omega_x)], [np.deg2rad(self.val_omega_y)], [np.deg2rad(self.val_omega_z)], [self.val_vel_x], [self.val_vel_y], [self.val_vel_z]])

            if np.linalg.norm(Twist_B[0:3]) < 1e-6:
                V = np.eye(3)
                e_VMatrix_theta = np.concatenate((np.concatenate((np.eye(3), V.dot(Twist_B[3:6]) * self.timer_value), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
            else:
                omegaSkew = np.array([[0, -Twist_B[2][0], Twist_B[1][0]], [Twist_B[2][0], 0, -Twist_B[0][0]], [-Twist_B[1][0], Twist_B[0][0], 0]]) * self.timer_value
                theta = np.linalg.norm(Twist_B[0:3] * self.timer_value)
                e_omegaSkew = np.eye(3) + ( np.sin(theta)/theta ) * omegaSkew + ( (1-np.cos(theta))/(theta**2) ) * np.linalg.matrix_power(omegaSkew, 2)
                V = np.eye(3) + ( (1-np.cos(theta))/(theta**2) ) * omegaSkew + ( (theta-np.sin(theta))/(theta**3) ) * np.linalg.matrix_power(omegaSkew, 2)
                e_VMatrix_theta = np.concatenate((np.concatenate((e_omegaSkew, V.dot(Twist_B[3:6] * self.timer_value)), axis=1), np.array([[0, 0, 0, 1]])), axis=0)
      
            B_T = e_VMatrix_theta.dot(B_T_zyx)  # Pre-multiply

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

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

        self.canvas.draw()

    def on_update_values(self):       
        self.B_p_x = float(self.ef_tx.text())
        self.B_p_y = float(self.ef_ty.text())
        self.B_p_z = float(self.ef_tz.text())
        self.B_p_qw = float(self.ef_qw.text())
        self.B_p_qx = float(self.ef_qx.text())
        self.B_p_qy = float(self.ef_qy.text())
        self.B_p_qz = float(self.ef_qz.text())

        if self.cbx_wx.isChecked():
            self.cbx_wx.setChecked(False)
            self.sld_wx.setValue(0.0)
        if self.cbx_wy.isChecked():
            self.cbx_wy.setChecked(False)
            self.sld_wy.setValue(0.0)
        if self.cbx_wz.isChecked():
            self.cbx_wz.setChecked(False)
            self.sld_wz.setValue(0.0) 
        if self.cbx_vx.isChecked():
            self.cbx_vx.setChecked(False)
            self.sld_vx.setValue(0.0)
        if self.cbx_vy.isChecked():
            self.cbx_vy.setChecked(False)
            self.sld_vy.setValue(0.0)
        if self.cbx_vz.isChecked():
            self.cbx_vz.setChecked(False)
            self.sld_vz.setValue(0.0)

        self.val_omega_x = self.sld_wx.value()
        self.val_omega_y = self.sld_wy.value()
        self.val_omega_z = self.sld_wz.value()
        self.val_vel_x = self.sld_vx.value()
        self.val_vel_y = self.sld_vy.value()
        self.val_vel_z = self.sld_vz.value()

        if self.rb_t0.isChecked() and self.scenario != 0:
            self.scenario = 0
            self.rb_t1.setChecked(False)
        elif self.rb_t1.isChecked() and self.scenario != 1:
            self.scenario = 1
            self.rb_t0.setChecked(False)

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
        
        self.cbx_wx = QCheckBox('reset')
        self.cbx_wx.setChecked(False)
        self.connect(self.cbx_wx, SIGNAL('stateChanged(int)'), self.on_update_values)
        
        self.sld_wx = DoubleSlider(Qt.Horizontal)
        self.sld_wx.setMinimum(-180.0)
        self.sld_wx.setMaximum(180.0)
        self.sld_wx.setValue(0.0)
        self.sld_wx.setTracking(True)
        self.sld_wx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_wx, SIGNAL('valueChanged(int)'), self.on_update_values)
        
        self.cbx_wy = QCheckBox('reset')
        self.cbx_wy.setChecked(False)
        self.connect(self.cbx_wy, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_wy = DoubleSlider(Qt.Horizontal)
        self.sld_wy.setMinimum(-180.0)
        self.sld_wy.setMaximum(180.0)
        self.sld_wy.setValue(0.0)
        self.sld_wy.setTracking(True)
        self.sld_wy.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_wy, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_wz = QCheckBox('reset')
        self.cbx_wz.setChecked(False)
        self.connect(self.cbx_wz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_wz = DoubleSlider(Qt.Horizontal)
        self.sld_wz.setMinimum(-180.0)
        self.sld_wz.setMaximum(180.0)
        self.sld_wz.setValue(0.0)
        self.sld_wz.setTracking(True)
        self.sld_wz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_wz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_vx = QCheckBox('reset')
        self.cbx_vx.setChecked(False)
        self.connect(self.cbx_vx, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_vx = DoubleSlider(Qt.Horizontal)
        self.sld_vx.setMinimum(-1.0)
        self.sld_vx.setMaximum(1.0)
        self.sld_vx.setValue(0.0)
        self.sld_vx.setTracking(True)
        self.sld_vx.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vx, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_vy = QCheckBox('reset')
        self.cbx_vy.setChecked(False)
        self.connect(self.cbx_vy, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_vy = DoubleSlider(Qt.Horizontal)
        self.sld_vy.setMinimum(-1.0)
        self.sld_vy.setMaximum(1.0)
        self.sld_vy.setValue(0.0)
        self.sld_vy.setTracking(True)
        self.sld_vy.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vy, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.cbx_vz = QCheckBox('reset')
        self.cbx_vz.setChecked(False)
        self.connect(self.cbx_vz, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.sld_vz = DoubleSlider(Qt.Horizontal)
        self.sld_vz.setMinimum(-1.0)
        self.sld_vz.setMaximum(1.0)
        self.sld_vz.setValue(0.0)
        self.sld_vz.setTracking(True)
        self.sld_vz.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vz, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.rb_t0 = QCheckBox('Twist #1')
        self.rb_t0.setChecked(True)
        self.connect(self.rb_t0, SIGNAL('stateChanged(int)'), self.on_update_values)

        self.rb_t1 = QCheckBox('Twist #2')
        self.rb_t1.setChecked(False)
        self.connect(self.rb_t1, SIGNAL('stateChanged(int)'), self.on_update_values)

        hbox_Rx = QHBoxLayout()
        for w in [ self.cbx_wx, QLabel('ωx'), QLabel('-180'), self.sld_wx, QLabel('180')]:
            hbox_Rx.addWidget(w)
            hbox_Rx.setAlignment(w, Qt.AlignVCenter)

        hbox_Ry = QHBoxLayout()
        for w in [ self.cbx_wy, QLabel('ωy'), QLabel('-180'), self.sld_wy, QLabel('180')]:
            hbox_Ry.addWidget(w)
            hbox_Ry.setAlignment(w, Qt.AlignVCenter)

        hbox_Rz = QHBoxLayout()
        for w in [ self.cbx_wz, QLabel('ωz'), QLabel('-180'), self.sld_wz, QLabel('180')]:
            hbox_Rz.addWidget(w)
            hbox_Rz.setAlignment(w, Qt.AlignVCenter)

        hbox_px = QHBoxLayout()
        for w in [ self.cbx_vx, QLabel('vx'), QLabel('-1'), self.sld_vx, QLabel('1')]:
            hbox_px.addWidget(w)
            hbox_px.setAlignment(w, Qt.AlignVCenter)

        hbox_py = QHBoxLayout()
        for w in [ self.cbx_vy, QLabel('vy'), QLabel('-1'), self.sld_vy, QLabel('1')]:
            hbox_py.addWidget(w)
            hbox_py.setAlignment(w, Qt.AlignVCenter)

        hbox_pz = QHBoxLayout()
        for w in [ self.cbx_vz, QLabel('vz'), QLabel('-1'), self.sld_vz, QLabel('1')]:
            hbox_pz.addWidget(w)
            hbox_pz.setAlignment(w, Qt.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ self.rb_t0, self.rb_t1, QLabel('t_x,y,z'), self.ef_tx, self.ef_ty, self.ef_tz, QLabel('q_w,x,y,z'), self.ef_qw, self.ef_qx, self.ef_qy, self.ef_qz]:
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

