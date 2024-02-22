import sys, os, random
import numpy as np

from transforms3d import *  #supersedes deprecated transformations.py
from transforms3d import _gohlketransforms


from PyQt6.QtCore import *
from PyQt6.QtGui import *


import matplotlib as plt

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from mpl_toolkits import mplot3d
from PyQt6.QtWidgets import *

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
        self.setWindowTitle('Rigid Body Motion: Transformations')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        #self.textbox.setText('1 2 3 4')


        self.val_rot_x = 0
        self.val_rot_y = 0
        self.val_rot_z = 0
        self.val_trans_x = 0
        self.val_trans_y = 0
        self.val_trans_z = 0

        self.scenario = 0

        self.cbx_px.setEnabled(False)
        self.sld_px.setEnabled(False)
        self.cbx_py.setEnabled(False)
        self.sld_py.setEnabled(False)
        self.cbx_pz.setEnabled(False)
        self.sld_pz.setEnabled(False)
        
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

        self.on_draw()

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_draw(self):
        self.axes.clear()        
        self.axes.grid(True)
        self.extraaxes.clear()

        if self.scenario == 0:
            # Rotation sequence around z-y-x-axes (post-multiply, each sequential rotation w.r.t. Body Frame)
            angle_z = np.deg2rad(self.val_rot_z)
            angle_y = np.deg2rad(self.val_rot_y)
            angle_x = np.deg2rad(self.val_rot_x)
            B_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0,                0,               1]])
            B_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            B_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])

            B_p_zyx = self.S_p
            B_R_zyx = ((self.S_R.dot(B_R_z)).dot(B_R_y)).dot(B_R_x)
            B_T_zyx = np.concatenate((np.concatenate((B_R_zyx, B_p_zyx), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            B_p = B_T_zyx[0:3,3]
            B_R = B_T_zyx[0:3,0:3]
            B_T = B_T_zyx

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_full = 1.0
            self.axes.set_title('Rotation', fontsize = 8)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-1.0, 1.0)
            self.axes.set_ylim(-1.0, 1.0)
            self.axes.set_zlim(-1.0, 1.0)

            B_q = _gohlketransforms.quaternion_from_matrix(B_R)
        
            self.extraaxes.set_title('Quaternion', fontsize = 8)
            self.extraaxes.set_xlim(-1.1, 1.1)
            self.extraaxes.barh(('q_w', 'q_x', 'q_y', 'q_z'), B_q, xerr=np.zeros(4), align='center')
            self.extraaxes.invert_yaxis()  # labels read top-to-bottom

            B_omega, B_theta = axangles.mat2axangle(B_R)
            if np.absolute(B_theta) > 1e-6:
                self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], 2.0*B_omega[0], 2.0*B_omega[1], 2.0*B_omega[2], color=['k'], arrow_length_ratio=0.0)


        elif self.scenario == 1:
            # Post-multiply, rotation w.r.t. Body Frame
            angle_z = np.deg2rad(45)
            angle_y = np.deg2rad(30)
            angle_x = np.deg2rad(15)
            B_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0, 0, 1]])
            B_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            B_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
        
            B_p_zyx = self.S_p
            B_R_zyx = ((self.S_R.dot(B_R_z)).dot(B_R_y)).dot(B_R_x)
            B_T_zyx = np.concatenate((np.concatenate((B_R_zyx, B_p_zyx), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            B_extra_R_z = axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(self.val_rot_z))
            B_extra_R_y = axangles.axangle2mat(np.array([0, 1, 0]), np.deg2rad(self.val_rot_y))
            B_extra_R_x = axangles.axangle2mat(np.array([1, 0, 0]), np.deg2rad(self.val_rot_x))

            B_p = B_p_zyx
            B_R = ((B_R_zyx.dot(B_extra_R_z)).dot(B_extra_R_y)).dot(B_extra_R_x)  # Pre-multiply
            B_T = np.concatenate((np.concatenate((B_R, B_p), axis=1), np.array([[0, 0, 0, 1]])), axis=0)

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.set_title('Composing w.r.t Body Frame: Post-Multiplying a z:[45deg],y:[30deg],x:[15deg] Rotation', fontsize = 8)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bx_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bx_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bx_zyx[2]-self.quiver_Bp_zyx[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_By_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_By_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_By_zyx[2]-self.quiver_Bp_zyx[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bz_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bz_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bz_zyx[2]-self.quiver_Bp_zyx[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-1.0, 1.0)
            self.axes.set_ylim(-1.0, 1.0)
            self.axes.set_zlim(-1.0, 1.0)

        elif self.scenario == 2:
            # Pre-multiply, rotation w.r.t. Space Frame
            angle_z = np.deg2rad(45)
            angle_y = np.deg2rad(30)
            angle_x = np.deg2rad(15)
            B_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0, 0, 1]])
            B_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            B_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
        
            B_p_zyx = self.S_p
            B_R_zyx = ((self.S_R.dot(B_R_z)).dot(B_R_y)).dot(B_R_x)
            B_T_zyx = np.concatenate((np.concatenate((B_R_zyx, B_p_zyx), axis=1), \
                                      np.array([[0, 0, 0, 1]])), axis=0)

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            B_extra_R_z = axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(self.val_rot_z))
            B_extra_R_y = axangles.axangle2mat(np.array([0, 1, 0]), np.deg2rad(self.val_rot_y))
            B_extra_R_x = axangles.axangle2mat(np.array([1, 0, 0]), np.deg2rad(self.val_rot_x))

            B_p = B_p_zyx
            B_R = B_extra_R_x.dot(B_extra_R_y.dot(B_extra_R_z.dot(B_R_zyx)))  # Pre-multiply
            B_T = np.concatenate((np.concatenate((B_R, B_p), axis=1), np.array([[0, 0, 0, 1]])), axis=0)

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.set_title('Composing w.r.t Space Frame: Pre-Multiplying a z:[45deg],y:[30deg],x:[15deg] Rotation', fontsize = 8)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bx_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bx_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bx_zyx[2]-self.quiver_Bp_zyx[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_By_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_By_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_By_zyx[2]-self.quiver_Bp_zyx[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bz_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bz_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bz_zyx[2]-self.quiver_Bp_zyx[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-1.0, 1.0)
            self.axes.set_ylim(-1.0, 1.0)
            self.axes.set_zlim(-1.0, 1.0)

        elif self.scenario == 3:
            # Translation w.r.t. all axes, Rotation sequence around z-y-x-axes (post-multiply, each sequential transformation w.r.t. Body Frame)
            B_p_zyx = np.array([[self.val_trans_x], \
                                [self.val_trans_y], \
                                [self.val_trans_z]])
            B_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), B_p_zyx), axis=1), \
                                            np.array([[0, 0, 0, 1]])), axis=0)
            angle_z = np.deg2rad(self.val_rot_z)
            angle_y = np.deg2rad(self.val_rot_y)
            angle_x = np.deg2rad(self.val_rot_x)
            B_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0,                0,               1]])
            B_T_z = np.concatenate((np.concatenate((B_R_z, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            B_T_y = np.concatenate((np.concatenate((B_R_y, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
            B_T_x = np.concatenate((np.concatenate((B_R_x, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)

            B_T_zyx = (((self.S_T.dot(B_T_zyx_trans)).dot(B_T_z)).dot(B_T_y)).dot(B_T_x)

            B_T = B_T_zyx
            B_p = B_T_zyx[0:3,3]
            B_R = B_T_zyx[0:3,0:3]

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_full = 1.0
            self.axes.set_title('Homogeneous Transformation', fontsize = 8)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-1.0, 1.0)
            self.axes.set_ylim(-1.0, 1.0)
            self.axes.set_zlim(-1.0, 1.0)

        elif self.scenario == 4:
            # Post-multiply, transformation w.r.t. Body Frame
            B_p_zyx = np.array([[0.50], \
                                [0.25], \
                                [0.75]])
            B_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), B_p_zyx), axis=1), \
                                            np.array([[0, 0, 0, 1]])), axis=0)
            angle_z = np.deg2rad(45)
            angle_y = np.deg2rad(30)
            angle_x = np.deg2rad(15)
            B_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0,                0,               1]])
            B_T_z = np.concatenate((np.concatenate((B_R_z, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            B_T_y = np.concatenate((np.concatenate((B_R_y, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
            B_T_x = np.concatenate((np.concatenate((B_R_x, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx = (((self.S_T.dot(B_T_zyx_trans)).dot(B_T_z)).dot(B_T_y)).dot(B_T_x)

            B_p_zyx = B_T_zyx[0:3,3]
            B_R_zyx = B_T_zyx[0:3,0:3]

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            B_extra_p = np.array([[self.val_trans_x], \
                                  [self.val_trans_y], \
                                  [self.val_trans_z]])
            B_extra_T_trans = np.concatenate((np.concatenate((np.eye(3), B_extra_p), axis=1), \
                                              np.array([[0, 0, 0, 1]])), axis=0)
            B_extra_R_z = axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(self.val_rot_z))
            B_extra_T_z = np.concatenate((np.concatenate((B_extra_R_z, np.zeros((3, 1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)
            B_extra_R_y = axangles.axangle2mat(np.array([0, 1, 0]), np.deg2rad(self.val_rot_y))
            B_extra_T_y = np.concatenate((np.concatenate((B_extra_R_y, np.zeros((3, 1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)
            B_extra_R_x = axangles.axangle2mat(np.array([1, 0, 0]), np.deg2rad(self.val_rot_x))
            B_extra_T_x = np.concatenate((np.concatenate((B_extra_R_x, np.zeros((3, 1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)

            B_T = (((B_T_zyx.dot(B_extra_T_trans)).dot(B_extra_T_z)).dot(B_extra_T_y)).dot(B_extra_T_x)  # Post-multiply

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.set_title('Composing w.r.t Body Frame: Post-Multiplying a Rot_zyx=[45,30,15]deg, Trans_xyz=[0.5,0.25,0.75]m', fontsize = 8)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bx_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bx_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bx_zyx[2]-self.quiver_Bp_zyx[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_By_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_By_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_By_zyx[2]-self.quiver_Bp_zyx[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bz_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bz_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bz_zyx[2]-self.quiver_Bp_zyx[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-1.0, 1.0)
            self.axes.set_ylim(-1.0, 1.0)
            self.axes.set_zlim(-1.0, 1.0)

        elif self.scenario == 5:
            # Pre-multiply, transformation w.r.t. Space Frame
            B_p_zyx = np.array([[0.50], \
                                [0.25], \
                                [0.75]])
            B_T_zyx_trans = np.concatenate((np.concatenate((np.eye(3), B_p_zyx), axis=1), \
                                            np.array([[0, 0, 0, 1]])), axis=0)
            angle_z = np.deg2rad(45)
            angle_y = np.deg2rad(30)
            angle_x = np.deg2rad(15)
            B_R_z = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], \
                              [np.sin(angle_z),  np.cos(angle_z), 0], \
                              [0,                0,               1]])
            B_T_z = np.concatenate((np.concatenate((B_R_z, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_R_y = np.array([[np.cos(angle_y),  0, np.sin(angle_y)], \
                              [0,                1,               0], \
                              [-np.sin(angle_y), 0, np.cos(angle_y)]])
            B_T_y = np.concatenate((np.concatenate((B_R_y, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_R_x = np.array([[1,               0,                0], \
                              [0, np.cos(angle_x), -np.sin(angle_x)], \
                              [0, np.sin(angle_x),  np.cos(angle_x)]])
            B_T_x = np.concatenate((np.concatenate((B_R_x, np.zeros((3, 1))), axis=1), \
                                    np.array([[0, 0, 0, 1]])), axis=0)
            B_T_zyx = (((self.S_T.dot(B_T_zyx_trans)).dot(B_T_z)).dot(B_T_y)).dot(B_T_x)

            B_p_zyx = B_T_zyx[0:3,3]
            B_R_zyx = B_T_zyx[0:3,0:3]

            self.quiver_Bp_zyx = B_T_zyx.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz_zyx = B_T_zyx.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            B_extra_p = np.array([[self.val_trans_x], \
                                  [self.val_trans_y], \
                                  [self.val_trans_z]])
            B_extra_T_trans = np.concatenate((np.concatenate((np.eye(3), B_extra_p), axis=1), \
                                              np.array([[0, 0, 0, 1]])), axis=0)
            B_extra_R_z = axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(self.val_rot_z))
            B_extra_T_z = np.concatenate((np.concatenate((B_extra_R_z, np.zeros((3, 1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)
            B_extra_R_y = axangles.axangle2mat(np.array([0, 1, 0]), np.deg2rad(self.val_rot_y))
            B_extra_T_y = np.concatenate((np.concatenate((B_extra_R_y, np.zeros((3, 1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)
            B_extra_R_x = axangles.axangle2mat(np.array([1, 0, 0]), np.deg2rad(self.val_rot_x))
            B_extra_T_x = np.concatenate((np.concatenate((B_extra_R_x, np.zeros((3, 1))), axis=1), \
                                          np.array([[0, 0, 0, 1]])), axis=0)

            B_T = B_extra_T_x.dot(B_extra_T_y.dot(B_extra_T_z.dot(B_extra_T_trans.dot(B_T_zyx))))  # Pre-multiply

            self.quiver_Bp = B_T.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = B_T.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = B_T.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = B_T.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_medium = 0.5
            scale_full = 1.0
            self.axes.set_title('Composing w.r.t Space Frame: Pre-Multiplying a Rot_zyx=[45,30,15]deg, Trans_xyz=[0.5,0.25,0.75]m', fontsize = 8)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], scale_small*self.quiver_Sx[2], color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], scale_small*self.quiver_Sy[2], color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.S_p[0], self.S_p[1], self.S_p[2], scale_small*self.quiver_Sz[0], scale_small*self.quiver_Sz[1], scale_small*self.quiver_Sz[2], color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bx_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bx_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bx_zyx[2]-self.quiver_Bp_zyx[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_By_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_By_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_By_zyx[2]-self.quiver_Bp_zyx[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp_zyx[0], self.quiver_Bp_zyx[1], self.quiver_Bp_zyx[2], scale_medium*(self.quiver_Bz_zyx[0]-self.quiver_Bp_zyx[0]), scale_medium*(self.quiver_Bz_zyx[1]-self.quiver_Bp_zyx[1]), scale_medium*(self.quiver_Bz_zyx[2]-self.quiver_Bp_zyx[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bx[2]-self.quiver_Bp[2]), color=['r'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), scale_full*(self.quiver_By[2]-self.quiver_Bp[2]), color=['g'], arrow_length_ratio=0.15)
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], self.quiver_Bp[2], scale_full*(self.quiver_Bz[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bz[1]-self.quiver_Bp[1]), scale_full*(self.quiver_Bz[2]-self.quiver_Bp[2]), color=['b'], arrow_length_ratio=0.15)
            self.axes.set_xlim(-1.0, 1.0)
            self.axes.set_ylim(-1.0, 1.0)
            self.axes.set_zlim(-1.0, 1.0)

        self.canvas.draw()

    def on_update_values(self):       
        if self.cbx_rx.isChecked():
            self.cbx_rx.setChecked(False)
            self.sld_rx.setValue(0.0)
        if self.cbx_ry.isChecked():
            self.cbx_ry.setChecked(False)
            self.sld_ry.setValue(0.0)
        if self.cbx_rz.isChecked():
            self.cbx_rz.setChecked(False)
            self.sld_rz.setValue(0.0) 
        if self.cbx_px.isChecked():
            self.cbx_px.setChecked(False)
            self.sld_px.setValue(0.0)
        if self.cbx_py.isChecked():
            self.cbx_py.setChecked(False)
            self.sld_py.setValue(0.0)
        if self.cbx_pz.isChecked():
            self.cbx_pz.setChecked(False)
            self.sld_pz.setValue(0.0)

        self.val_rot_x = self.sld_rx.value()
        self.val_rot_y = self.sld_ry.value()
        self.val_rot_z = self.sld_rz.value()
        self.val_trans_x = self.sld_px.value()
        self.val_trans_y = self.sld_py.value()
        self.val_trans_z = self.sld_pz.value()

        if self.rb_r0.isChecked() and self.scenario != 0:
            self.scenario = 0
            self.rb_r1.setChecked(False)
            self.rb_r2.setChecked(False)
            self.rb_p0.setChecked(False)
            self.rb_p1.setChecked(False)
            self.rb_p2.setChecked(False)

            self.cbx_px.setEnabled(False)
            self.sld_px.setEnabled(False)
            self.cbx_py.setEnabled(False)
            self.sld_py.setEnabled(False)
            self.cbx_pz.setEnabled(False)
            self.sld_pz.setEnabled(False)

            self.extraaxes.set_visible(True)
        elif self.rb_r1.isChecked() and self.scenario != 1:
            self.scenario = 1
            self.rb_r0.setChecked(False)
            self.rb_r2.setChecked(False)
            self.rb_p0.setChecked(False)
            self.rb_p1.setChecked(False)
            self.rb_p2.setChecked(False)

            self.cbx_px.setEnabled(False)
            self.sld_px.setEnabled(False)
            self.cbx_py.setEnabled(False)
            self.sld_py.setEnabled(False)
            self.cbx_pz.setEnabled(False)
            self.sld_pz.setEnabled(False)

            self.extraaxes.set_visible(False)
        elif self.rb_r2.isChecked() and self.scenario != 2:
            self.scenario = 2
            self.rb_r0.setChecked(False)
            self.rb_r1.setChecked(False)
            self.rb_p0.setChecked(False)
            self.rb_p1.setChecked(False)
            self.rb_p2.setChecked(False)

            self.cbx_px.setEnabled(False)
            self.sld_px.setEnabled(False)
            self.cbx_py.setEnabled(False)
            self.sld_py.setEnabled(False)
            self.cbx_pz.setEnabled(False)
            self.sld_pz.setEnabled(False)

            self.extraaxes.set_visible(False)
        elif self.rb_p0.isChecked() and self.scenario != 3:
            self.scenario = 3
            self.rb_r0.setChecked(False)
            self.rb_r1.setChecked(False)
            self.rb_r2.setChecked(False)
            self.rb_p1.setChecked(False)
            self.rb_p2.setChecked(False)

            self.cbx_px.setEnabled(True)
            self.sld_px.setEnabled(True)
            self.cbx_py.setEnabled(True)
            self.sld_py.setEnabled(True)
            self.cbx_pz.setEnabled(True)
            self.sld_pz.setEnabled(True)

            self.extraaxes.set_visible(False)
        elif self.rb_p1.isChecked() and self.scenario != 4:
            self.scenario = 4
            self.rb_r0.setChecked(False)
            self.rb_r1.setChecked(False)
            self.rb_r2.setChecked(False)
            self.rb_p0.setChecked(False)
            self.rb_p2.setChecked(False)

            self.cbx_px.setEnabled(True)
            self.sld_px.setEnabled(True)
            self.cbx_py.setEnabled(True)
            self.sld_py.setEnabled(True)
            self.cbx_pz.setEnabled(True)
            self.sld_pz.setEnabled(True)

            self.extraaxes.set_visible(False)
        elif self.rb_p2.isChecked() and self.scenario != 5:
            self.scenario = 5
            self.rb_r0.setChecked(False)
            self.rb_r1.setChecked(False)
            self.rb_r2.setChecked(False)
            self.rb_p0.setChecked(False)
            self.rb_p1.setChecked(False)

            self.cbx_px.setEnabled(True)
            self.sld_px.setEnabled(True)
            self.cbx_py.setEnabled(True)
            self.sld_py.setEnabled(True)
            self.cbx_pz.setEnabled(True)
            self.sld_pz.setEnabled(True)

            self.extraaxes.set_visible(False)

        self.on_draw()

    def create_main_frame(self):
        self.main_frame = QWidget()

        self.dpi = 100
        self.fig = Figure((5.0, 15.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        self.axes = self.fig.add_subplot(111, projection='3d', proj_type='ortho')
        self.axes.set_aspect('equal')
        self.extraaxes = self.fig.add_subplot(241)
        self.extraaxes.set_visible(True)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        

        #self.textbox = QLineEdit()
        #self.textbox.setMinimumWidth(200)
        #self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_update_values)
        
        #self.draw_button = QPushButton("&Draw")
        #self.connect(self.draw_button, SIGNAL('clicked()'), self.on_update_values)
        
        self.cbx_rx = QCheckBox('reset')
        self.cbx_rx.setChecked(False)
        self.cbx_rx.stateChanged.connect(self.on_update_values)
        
        self.sld_rx = DoubleSlider(Qt.Orientation.Horizontal)
        self.sld_rx.setMinimum(-180.0)
        self.sld_rx.setMaximum(180.0)
        self.sld_rx.setValue(0.0)
        self.sld_rx.setTracking(True)
        self.sld_rx.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.sld_rx.valueChanged.connect(self.on_update_values)
        
        self.cbx_ry = QCheckBox('reset')
        self.cbx_ry.setChecked(False)
        self.cbx_ry.stateChanged.connect(self.on_update_values)

        self.sld_ry = DoubleSlider(Qt.Orientation.Horizontal)
        self.sld_ry.setMinimum(-180.0)
        self.sld_ry.setMaximum(180.0)
        self.sld_ry.setValue(0.0)
        self.sld_ry.setTracking(True)
        self.sld_ry.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.sld_ry.valueChanged.connect(self.on_update_values)

        self.cbx_rz = QCheckBox('reset')
        self.cbx_rz.setChecked(False)
        self.cbx_rz.stateChanged.connect(self.on_update_values)

        self.sld_rz = DoubleSlider(Qt.Orientation.Horizontal)
        self.sld_rz.setMinimum(-180.0)
        self.sld_rz.setMaximum(180.0)
        self.sld_rz.setValue(0.0)
        self.sld_rz.setTracking(True)
        self.sld_rz.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.sld_rz.valueChanged.connect(self.on_update_values)

        self.cbx_px = QCheckBox('reset')
        self.cbx_px.setChecked(False)
        self.cbx_px.stateChanged.connect(self.on_update_values)

        self.sld_px = DoubleSlider(Qt.Orientation.Horizontal)
        self.sld_px.setMinimum(-1.0)
        self.sld_px.setMaximum(1.0)
        self.sld_px.setValue(0.0)
        self.sld_px.setTracking(True)
        self.sld_px.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.sld_px.valueChanged.connect(self.on_update_values)

        self.cbx_py = QCheckBox('reset')
        self.cbx_py.setChecked(False)
        self.cbx_py.stateChanged.connect(self.on_update_values)

        self.sld_py = DoubleSlider(Qt.Orientation.Horizontal)
        self.sld_py.setMinimum(-1.0)
        self.sld_py.setMaximum(1.0)
        self.sld_py.setValue(0.0)
        self.sld_py.setTracking(True)
        self.sld_py.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.sld_py.valueChanged.connect(self.on_update_values)

        self.cbx_pz = QCheckBox('reset')
        self.cbx_pz.setChecked(False)
        self.cbx_pz.stateChanged.connect(self.on_update_values) 

        self.sld_pz = DoubleSlider(Qt.Orientation.Horizontal)
        self.sld_pz.setMinimum(-1.0)
        self.sld_pz.setMaximum(1.0)
        self.sld_pz.setValue(0.0)
        self.sld_pz.setTracking(True)
        self.sld_pz.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.sld_pz.valueChanged.connect(self.on_update_values)

        self.rb_r0 = QCheckBox('Rot #1')
        self.rb_r0.setChecked(True)
        self.rb_r0.stateChanged.connect(self.on_update_values)

        self.rb_r1 = QCheckBox('Rot #2')
        self.rb_r1.setChecked(False)
        self.rb_r1.stateChanged.connect(self.on_update_values)

        self.rb_r2 = QCheckBox('Rot #3')
        self.rb_r2.setChecked(False)
        self.rb_r2.stateChanged.connect(self.on_update_values)

        self.rb_p0 = QCheckBox('Homog #1')
        self.rb_p0.setChecked(False)
        self.rb_p0.stateChanged.connect(self.on_update_values)

        self.rb_p1 = QCheckBox('Homog #2')
        self.rb_p1.setChecked(False)
        self.rb_p1.stateChanged.connect(self.on_update_values)

        self.rb_p2 = QCheckBox('Homog #3')
        self.rb_p2.setChecked(False)
        self.rb_p2.stateChanged.connect(self.on_update_values)

        hbox_Rx = QHBoxLayout()
        for w in [ self.cbx_rx, QLabel('Rx'), QLabel('-180'), self.sld_rx, QLabel('180')]:
            hbox_Rx.addWidget(w)
            hbox_Rx.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

        hbox_Ry = QHBoxLayout()
        for w in [ self.cbx_ry, QLabel('Ry'), QLabel('-180'), self.sld_ry, QLabel('180')]:
            hbox_Ry.addWidget(w)
            hbox_Ry.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

        hbox_Rz = QHBoxLayout()
        for w in [ self.cbx_rz, QLabel('Rz'), QLabel('-180'), self.sld_rz, QLabel('180')]:
            hbox_Rz.addWidget(w)
            hbox_Rz.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

        hbox_px = QHBoxLayout()
        for w in [ self.cbx_px, QLabel('px'), QLabel('-1'), self.sld_px, QLabel('1')]:
            hbox_px.addWidget(w)
            hbox_px.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

        hbox_py = QHBoxLayout()
        for w in [ self.cbx_py, QLabel('py'), QLabel('-1'), self.sld_py, QLabel('1')]:
            hbox_py.addWidget(w)
            hbox_py.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

        hbox_pz = QHBoxLayout()
        for w in [ self.cbx_pz, QLabel('pz'), QLabel('-1'), self.sld_pz, QLabel('1')]:
            hbox_pz.addWidget(w)
            hbox_pz.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ self.rb_r0, self.rb_r1, self.rb_r2, self.rb_p0, self.rb_p1, self.rb_p2]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignmentFlag.AlignVCenter)

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
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec()

if __name__ == "__main__":
    main()
