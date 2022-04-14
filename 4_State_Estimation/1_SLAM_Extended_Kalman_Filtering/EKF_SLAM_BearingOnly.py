import sys, os, random
import numpy as np
import time

from transforms3d import *  #supersedes deprecated transformations.py

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import cm

from mpl_toolkits import mplot3d

from itertools import product, combinations

from numpy.random import default_rng


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


def randomSign():
    return 1 if random.random() < 0.5 else -1

def wrapToPi(radians):
    while radians >= np.pi:
        radians = radians - 2.0*np.pi
    while radians < -np.pi:
        radians = radians + 2.0*np.pi
    return radians

def setupBody(axles_distance, half_length, half_width, half_height):
    body_vertices = list()
    body_vertices.append([[-half_length+axles_distance/2, -half_width, -half_height], [-half_length+axles_distance/2, -half_width, half_height]])
    body_vertices.append([[-half_length+axles_distance/2, -half_width, -half_height], [half_length+axles_distance/2, -half_width, -half_height]])
    body_vertices.append([[-half_length+axles_distance/2, -half_width, -half_height], [half_length+axles_distance/2, -half_width, -half_height]])
    body_vertices.append([[-half_length+axles_distance/2, -half_width, half_height], [-half_length+axles_distance/2, half_width, half_height]])
    body_vertices.append([[-half_length+axles_distance/2, -half_width, half_height], [half_length+axles_distance/2, -half_width, half_height]])
    body_vertices.append([[-half_length+axles_distance/2, half_width, -half_height], [-half_length+axles_distance/2, half_width, half_height]])
    body_vertices.append([[-half_length+axles_distance/2, half_width, -half_height], [half_length+axles_distance/2, half_width, -half_height]])
    body_vertices.append([[-half_length+axles_distance/2, half_width, half_height], [half_length+axles_distance/2, half_width, half_height]])
    body_vertices.append([[half_length+axles_distance/2, -half_width, -half_height], [half_length+axles_distance/2, -half_width, half_height]])
    body_vertices.append([[half_length+axles_distance/2, -half_width, -half_height], [half_length+axles_distance/2, half_width, -half_height]])
    body_vertices.append([[half_length+axles_distance/2, -half_width, half_height], [half_length+axles_distance/2, half_width, half_height]])
    body_vertices.append([[half_length+axles_distance/2, half_width, -half_height], [half_length+axles_distance/2, half_width, half_height]])
    return body_vertices

def colormap(index):
    if index > 0:
        selected_color = "m"
    elif index == 0:
        selected_color = "k"
    elif index == -1:
        selected_color = "b"
    elif index == -2:
        selected_color = "y"
    elif index == -3:
        selected_color = "r" 
    return selected_color

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Non-Holonomic Robot')

        self.m_groundtruth = np.transpose(np.array([[3.0, 2.0], \
                                                   [5.0, -1.0], \
                                                   [-5.0, -1.0], \
                                                   [-4.0, -5.0], \
                                                   [-3.0, 5.0]]))
        self.num_landmarks = len(np.transpose(self.m_groundtruth))
        self.landmarks_quality_history = [-1] * self.num_landmarks  # negative: uninintialized , 0-threshold: valid , >threshold: replace
        self.landmarks_quality_counter_threshold = 5

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        qApp.installEventFilter(self)

        self.fire_event = False

        self.rng = default_rng()

        for i in range(self.num_landmarks): 
            self.m_groundtruth_text[i].setText(str(self.m_groundtruth[0][i])+","+str(self.m_groundtruth[1][i]))

        ksi_0 = np.array([[0.0], \
                          [0.0], \
                          [np.deg2rad(0.0)]])  # [x, y, theta]
        self.ksi_groundtruth = np.copy(ksi_0)
        self.ksi_estim = np.copy(ksi_0)
        self.ksi_hat_previous = np.copy(ksi_0)

        self.m_estim = np.copy(self.m_groundtruth)
        for i in range(self.num_landmarks):
            #self.m_estim[:,[i]] = 5.0 * self.m_groundtruth[:,[i]] / np.linalg.norm(self.m_groundtruth[:,[i]])  #initialize at 5m range (blindly)
            z_mi_dxy_init = self.m_groundtruth[:,[i]] - self.ksi_groundtruth[0:1 + 1]
            z_mi_dr_init = np.linalg.norm(z_mi_dxy_init) + self.rng.normal(0.0, 1.0)  #initialize at a depth of +-1.0m sigma (assumes we have some disparity-based delayed initialization pipeline for depth)
            z_mi_dtheta_init = np.arctan2(z_mi_dxy_init[1][0], z_mi_dxy_init[0][0]) - self.ksi_groundtruth[2][0] + self.rng.normal(0.0,np.deg2rad(5.0))  #initialize at a bearing of +-5.0deg sigma
            z_mi_dtheta_init = wrapToPi(z_mi_dtheta_init)

            self.m_estim[:,[i]] = np.array([[self.ksi_estim[0][0] + z_mi_dr_init*np.cos(self.ksi_estim[2][0] + z_mi_dtheta_init)], \
                                            [self.ksi_estim[1][0] + z_mi_dr_init*np.sin(self.ksi_estim[2][0] + z_mi_dtheta_init)]])

        self.u_t = np.array([[0.0], \
                             [0.0]])  # input forward and turn velocity

        self.z_mi = np.empty((2,self.num_landmarks))  # measurements vector for all landmarks
        self.z_mi[:,:] = np.NaN

        for i in range(self.num_landmarks):
            self.landmarks_quality_history[i] = -1

        self.timer_period = 0.10
        self.timer_value = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_timer)
        self.timer.start(self.timer_period * 1000)

        #r = [-1, 1]
        #for s, e in combinations(np.array(list(product(r, r, r))), 2):
        #    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        #        print(s, e)
        #        self.axes.plot3D(*zip(s, e), color="k")
        self.axles_distance = 0
        self.half_length = 0.5
        self.half_width = 0.5
        self.half_height = 0.10
        self.body_vertices = setupBody(self.axles_distance, self.half_length, self.half_width, self.half_height)

        self.sigma_Rwheel = 0.10  # wheel velocity
        self.sigma_Lwheel = 0.10  # wheel velocity

        self.sigma_v = np.sqrt( ((0.5)**2) * ((self.sigma_Rwheel)**2) + ((0.5)**2) * ((self.sigma_Lwheel)**2) )  # forward velocity
        self.sigma_w = np.sqrt( ((1.0/(2*self.half_width))**2) * ((self.sigma_Rwheel)**2) + ((-1.0/(2*self.half_width))**2) * ((self.sigma_Lwheel)**2) )  # rotational velocity

        self.sigma_bearing = np.deg2rad(5.0)  # landmark relative bearing # static,if changing then for every pose-landmark estimate based on range
   
        self.Q = 1.0 * np.array([[self.sigma_v**2,               0], \
                                 [0,               self.sigma_w**2]])
        self.R = 1.0 * np.array([[self.sigma_bearing**2]])

        S_ksi_estim = 1.0 * np.array([[0.00**2,       0,                  0], \
                                      [0,       0.00**2,                  0], \
                                      [0,             0, np.deg2rad(0.0)**2]])
        self.W_S_ksi_trans_init, self.V_S_ksi_trans_init = np.linalg.eig(S_ksi_estim[0 : 1 + 1 , 0 : 1 +1])
        idx_S_ksi_trans_init = self.W_S_ksi_trans_init.argsort()[::-1]   
        self.W_S_ksi_trans_init = self.W_S_ksi_trans_init[idx_S_ksi_trans_init]
        self.V_S_ksi_trans_init = self.V_S_ksi_trans_init[:,idx_S_ksi_trans_init]
        self.W_S_ksi_rot_init, self.V_S_ksi_rot_init = np.linalg.eig(S_ksi_estim[2 : 2 + 1 , 2 : 2 +1])
        idx_S_ksi_rot_init = self.W_S_ksi_rot_init.argsort()[::-1]   
        self.W_S_ksi_rot_init = self.W_S_ksi_rot_init[idx_S_ksi_rot_init]
        self.V_S_ksi_rot_init = self.V_S_ksi_rot_init[:,idx_S_ksi_rot_init]

        self.S_m_new = 1.0 * np.array([[10000.0**2,          0], \
                                       [0,          10000.0**2]])
        self.S_estim = np.concatenate((np.concatenate((S_ksi_estim, np.zeros((3, self.m_groundtruth.size))), axis=1), \
                                       np.concatenate((np.transpose(np.zeros((3, self.m_groundtruth.size))), np.kron(np.eye(self.num_landmarks,dtype=float),self.S_m_new)), axis=1)), axis=0)

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

        #self.on_draw()

    def eventFilter(self, source, event):
        if event.type() == QEvent.KeyPress:
            #print('KeyPress: %s [%r]' % (event.key(), source))
            if event.key() == Qt.Key_Up:  # up
                self.u_t[0][0] = self.u_t[0][0] + 0.25
            elif event.key() == Qt.Key_Down:  # down
                self.u_t[0][0] = self.u_t[0][0] - 0.25
            elif event.key() == Qt.Key_Left:  # left
                self.u_t[1][0] = self.u_t[1][0] + 0.1
            elif event.key() == Qt.Key_Right:  # right
                self.u_t[1][0] = self.u_t[1][0] - 0.1
            elif event.key() == Qt.Key_Space:  # spacebar
                self.fire_event = not self.fire_event

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

        if np.abs(self.u_t[0][0]) > 0:
            if np.abs(self.u_t[0][0]) <= 0.05:
                self.u_t[0][0] = 0
            else:
                self.u_t[0][0] = self.u_t[0][0] - np.sign(self.u_t[0][0]) * 0.05
        if np.abs(self.u_t[1][0]) > 0:
            if np.abs(self.u_t[1][0]) <= 0.01:
                self.u_t[1][0] = 0
            else:
                self.u_t[1][0] = self.u_t[1][0] - np.sign(self.u_t[1][0]) * 0.01

        # update low-level Joint-Space commands (wheels) & add noise
        u_R = self.u_t[0][0] + self.u_t[1][0] * self.half_width + np.max(np.array([-2.0*self.sigma_Rwheel, np.min(np.array([2.0*self.sigma_Rwheel, self.rng.normal(0.0, self.sigma_Rwheel)]))]))
        u_L = self.u_t[0][0] - self.u_t[1][0] * self.half_width + np.max(np.array([-2.0*self.sigma_Lwheel, np.min(np.array([2.0*self.sigma_Lwheel, self.rng.normal(0.0, self.sigma_Lwheel)]))]))

        # calculate Body Frame kinematics based on joint-space velocities
        u_vel = (0.5)*(u_R + u_L);
        u_rot = (1.0/(2*self.half_width))*(u_R - u_L)

        x_dot = u_vel * np.cos( self.ksi_groundtruth[2][0] )
        y_dot = u_vel * np.sin( self.ksi_groundtruth[2][0] )
        theta_dot = u_rot;  

        # update Space Frame kinematics 
        self.ksi_groundtruth[0][0] = self.ksi_groundtruth[0][0] + x_dot * self.timer_period
        self.ksi_groundtruth[1][0] = self.ksi_groundtruth[1][0] + y_dot * self.timer_period
        self.ksi_groundtruth[2][0] = self.ksi_groundtruth[2][0] + theta_dot * self.timer_period
        #wrapToPi for theta_groundtruth
        #self.ksi_groundtruth[2][0] = wrapToPi(self.ksi_groundtruth[2][0])

        #print(self.timer_value)

        self.on_draw()

    def on_draw(self):
        self.axes.clear()        
        self.axes.grid(True)

        # draw groundtruth
        T_SB = np.concatenate((np.concatenate((axangles.axangle2mat(np.array([0, 0, 1]), self.ksi_groundtruth[2][0]), np.array([[self.ksi_groundtruth[0][0]], [self.ksi_groundtruth[1][0]], [0]])), axis=1),
                               np.array([[0, 0, 0, 1]])), axis=0)

        self.quiver_Bp = T_SB.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
        self.quiver_Bx = T_SB.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
        self.quiver_By = T_SB.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
        self.quiver_Bz = T_SB.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

        scale = 15.0
        self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], 1.0*(self.quiver_Bx[0]-self.quiver_Bp[0]), 1.0*(self.quiver_Bx[1]-self.quiver_Bp[1]), color=['r'], scale=scale, width=0.005)
        self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], 1.0*(self.quiver_By[0]-self.quiver_Bp[0]), 1.0*(self.quiver_By[1]-self.quiver_Bp[1]), color=['g'], scale=scale, width=0.005)
        self.axes.set_xlim(self.quiver_Bp[0]-scale, self.quiver_Bp[0]+scale)
        self.axes.set_ylim(self.quiver_Bp[1]-scale, self.quiver_Bp[1]+scale)
        #self.axes.set_xlim(-scale, scale)
        #self.axes.set_ylim(-scale, scale)

        # draw groundtruth body
        for vertex_pair in self.body_vertices:
            v0 = vertex_pair[0] 
            v1 = vertex_pair[1]
            v0_transformed = T_SB.dot(np.array([[v0[0]], [v0[1]], [v0[2]], [1]]))
            v1_transformed = T_SB.dot(np.array([[v1[0]], [v1[1]], [v1[2]], [1]]))
            #self.axes.plot3D(*zip(v0_transformed[0:2], v1_transformed[0:2]), color="b")
            self.axes.plot([v0_transformed[0], v1_transformed[0]], [v0_transformed[1], v1_transformed[1]], color="b")

        # form range bearing measurements, corrupted by noise
        z_mi = np.empty((1,self.num_landmarks)) 
        z_mi[:,:] = np.NaN
        for i in range(self.num_landmarks):   
            z_mi_dxy_groundtruth = self.m_groundtruth[:,[i]] - self.ksi_groundtruth[0:1 + 1]
            z_mi_dtheta_groundtruth = np.arctan2(z_mi_dxy_groundtruth[1][0], z_mi_dxy_groundtruth[0][0]) - self.ksi_groundtruth[2][0] + np.max(np.array([-2.0*self.sigma_bearing, np.min(np.array([2.0*self.sigma_bearing, self.rng.normal(0.0, self.sigma_bearing)]))]))
            #wrapToPi for bearing measurement
            z_mi_dtheta_groundtruth = wrapToPi(z_mi_dtheta_groundtruth)
            z_mi[:,[i]] = np.array([[z_mi_dtheta_groundtruth]]) 

        # draw groundtruth landmarks
        for i in range(self.num_landmarks):
            self.axes.scatter([self.m_groundtruth[0][i]], [self.m_groundtruth[1][i]], color="k")
            #self.axes.plot([self.ksi_groundtruth[0][0], self.m_groundtruth[0][i]], [self.ksi_groundtruth[1][0], self.m_groundtruth[1][i]], color="k", linestyle="dashed", linewidth=1.0)
            self.axes.plot([self.ksi_groundtruth[0][0], self.ksi_groundtruth[0][0] + 100.0*np.cos(z_mi[0][i]+self.ksi_groundtruth[2][0])], [self.ksi_groundtruth[1][0], self.ksi_groundtruth[1][0] + 100.0*np.sin(z_mi[0][i]+self.ksi_groundtruth[2][0])], color="k", linestyle="dashed", linewidth=1.0)


        ### Extended Kalman Filtering ###
        #if self.fire_event == False:
        #    self.canvas.draw()
        #    return
        print("EKF begin")

        # predicted estimate
        ksi_hat = np.copy(self.ksi_estim)
        ksi_hat[0][0] = ksi_hat[0][0] + self.timer_period * self.u_t[0][0]*np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)
        ksi_hat[1][0] = ksi_hat[1][0] + self.timer_period * self.u_t[0][0]*np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)
        ksi_hat[2][0] = ksi_hat[2][0] + self.timer_period * self.u_t[1][0]

        F_ksi = np.array([[1, 0, self.timer_period * -self.u_t[0][0]*np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                          [0, 1, self.timer_period *  self.u_t[0][0]*np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                          [0, 0, 1]])

        F_u = np.array([[self.timer_period * np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period), -(0.5*self.timer_period) * self.timer_period * self.u_t[0][0]*np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                        [self.timer_period * np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period),  (0.5*self.timer_period) * self.timer_period * self.u_t[0][0]*np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                        [0                                                                                            ,  self.timer_period]])

        Fksi_Sksiksi_FksiT = F_ksi.dot(self.S_estim[0 : 2 + 1, 0 : 2 + 1]).dot(np.transpose(F_ksi)) + F_u.dot(self.Q).dot(np.transpose(F_u))
        Fksi_Sksim = F_ksi.dot(self.S_estim[0 : 2 + 1, 3 + 2*0 : 3+ 2*self.num_landmarks+1 + 1])

        # update pose estimate
        self.ksi_estim = np.copy(ksi_hat)

        # update pose - landmarks covariance
        self.S_estim[0 : 2 + 1                               , 0 : 2 + 1]                               = Fksi_Sksiksi_FksiT
        self.S_estim[0 : 2 + 1                               , 3 + 2*0 : 3+ 2*self.num_landmarks+1 + 1] = Fksi_Sksim
        self.S_estim[3 + 2*0 : 3+ 2*self.num_landmarks+1 + 1 , 0 : 2 + 1]                               = np.transpose(Fksi_Sksim)

        # # ARTIFICIAL COVARIANCE INFLATION
        # W_trans, V_trans = np.linalg.eig(self.S_estim[0 : 1 + 1 , 0 : 1 +1])
        # idx_trans = W_trans.argsort()[::-1]   
        # W_trans = W_trans[idx_trans]
        # V_trans = V_trans[:,idx_trans]
        # if W_trans[1]<self.W_S_ksi_trans_init[1]:
        #     W_trans_inflated = np.copy(W_trans) 
        #     W_trans_inflated[1] = self.W_S_ksi_trans_init[1]
        #     if W_trans[0]<self.W_S_ksi_trans_init[0]:
        #         W_trans_inflated[0] = self.W_S_ksi_trans_init[0]
        #     D_trans = np.diag(W_trans_inflated)    
        #     self.S_estim[0 : 1 + 1, 0 : 1 + 1] = V_trans.dot(D_trans).dot(np.linalg.pinv(V_trans))

        landmarks_quality = [0] * self.num_landmarks
        # form range bearing predictions and Jacobian entries, conduct Landmark i update
        for i in range(self.num_landmarks):  
            # Error
            z_mi_dxy_hat = self.m_estim[:,[i]] - self.ksi_estim[0 : 1 + 1,:] 
            z_mi_dr_hat = np.linalg.norm(z_mi_dxy_hat)
            z_mi_dtheta_hat = np.arctan2(z_mi_dxy_hat[1][0], z_mi_dxy_hat[0][0]) - self.ksi_estim[2][0]
            #wrapToPi for bearing delta (we don't do actual SE(2) math)
            z_mi_dtheta_hat = wrapToPi(z_mi_dtheta_hat)
            z_mi_hat = np.array([[z_mi_dtheta_hat]]) 

            # Check distance
            if z_mi_dr_hat > 15.0:
                print(i, " dist")
                landmarks_quality[i] = -1
                if self.landmarks_quality_history[i] >= 0:
                    self.landmarks_quality_history[i] = self.landmarks_quality_history[i] + 1
                continue

            # Check inverse distance to avoid an ill-conditioned Jacobian
            if z_mi_dr_hat < 0.1:
                print(i, " range")
                landmarks_quality[i] = -1
                continue

            z_mi_dx_hat = z_mi_dxy_hat[0][0]
            z_mi_dy_hat = z_mi_dxy_hat[1][0]
            q = 1.0/z_mi_dr_hat
            q_sq = q**2 
            H_mi = np.concatenate((np.concatenate((np.concatenate((np.array([[z_mi_dy_hat*q_sq, -z_mi_dx_hat*q_sq, -1]]), \
                                                                                                                          np.zeros((1,2*i))),axis=1), np.array([[-z_mi_dy_hat*q_sq, z_mi_dx_hat*q_sq]])), axis=1), np.zeros((1,2*(-1+self.num_landmarks-i)))), axis=1)

            H_Shat_Ht = H_mi.dot(self.S_estim).dot(np.transpose(H_mi))

            R_mi = H_Shat_Ht + self.R

            R_mi_inv = np.linalg.pinv(R_mi, hermitian=True) 

            z_err = z_mi[:,[i]] - z_mi_hat
            #wrapToPi for bearing error (we don't do actual SE(2) math)
            z_err[0][0] = wrapToPi(z_err[0][0])

            # Check Mahalanobis distance
            Mah_i = np.transpose(z_err).dot(R_mi_inv.dot(z_err))[0][0]
            if Mah_i > 9.0 or Mah_i < 0.0:
                print(i, " Mah: ", Mah_i)
                landmarks_quality[i] = -2
                if self.landmarks_quality_history[i] >= 0:
                    self.landmarks_quality_history[i] = self.landmarks_quality_history[i] + 1
                continue

            x_hat = np.copy(self.ksi_estim)
            for j in range(self.num_landmarks):
                x_hat = np.concatenate((x_hat,
                                        self.m_estim[:,[j]]), axis=0)

            # EKF Update step
            K_mi = self.S_estim.dot(np.transpose(H_mi)).dot(R_mi_inv)
            x_upd = x_hat + K_mi.dot(z_err)

            S_upd = ( np.eye(3+2*self.num_landmarks) - K_mi.dot(H_mi) ).dot(self.S_estim);
            #S_upd = S_hat - K_mi.dot(R_mi).dot(np.transpose(K_mi))

            # Check Positive-Definiteness
            if not np.all(np.linalg.eigvals(S_upd) > 0):
                print(i, " eig: "), print(S_upd)
                landmarks_quality[i] = -3
                if self.landmarks_quality_history[i] >= 0:
                    self.landmarks_quality_history[i] = self.landmarks_quality_history[i] + 1
                continue

            # Update pose estimate
            self.ksi_estim = x_upd[0 : 2 + 1]

            # Update landmark estimate
            for j in range(self.num_landmarks):
                self.m_estim[:,[j]] = x_upd[3 + 2*j : 3 + 2*j+1 + 1]

            # Update pose - landmarks covariance
            self.S_estim = np.copy(S_upd)

            print(i, " ok")
            landmarks_quality[i] = 1
            self.landmarks_quality_history[i] = 0

        # draw estimate
        T_SB_estim = np.concatenate((np.concatenate((axangles.axangle2mat(np.array([0, 0, 1]), self.ksi_estim[2][0]), np.array([[self.ksi_estim[0][0]], [self.ksi_estim[1][0]], [0]])), axis=1),
                                     np.array([[0, 0, 0, 1]])), axis=0)

        self.quiver_Bp_estim = T_SB_estim.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
        self.quiver_Bx_estim = T_SB_estim.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
        self.quiver_By_estim = T_SB_estim.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
        self.quiver_Bz_estim = T_SB_estim.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

        scale = 15.0
        self.axes.quiver(self.quiver_Bp_estim[0], self.quiver_Bp_estim[1], 1.0*(self.quiver_Bx_estim[0]-self.quiver_Bp_estim[0]), 1.0*(self.quiver_Bx_estim[1]-self.quiver_Bp_estim[1]), color=['r'], scale=scale, width=0.005)
        self.axes.quiver(self.quiver_Bp_estim[0], self.quiver_Bp_estim[1], 1.0*(self.quiver_By_estim[0]-self.quiver_Bp_estim[0]), 1.0*(self.quiver_By_estim[1]-self.quiver_Bp_estim[1]), color=['g'], scale=scale, width=0.005)
        self.axes.set_xlim(self.quiver_Bp[0]-scale, self.quiver_Bp[0]+scale)
        self.axes.set_ylim(self.quiver_Bp[1]-scale, self.quiver_Bp[1]+scale)
        #self.axes.set_xlim(-scale, scale)
        #self.axes.set_ylim(-scale, scale)

        # draw estimate body
        for vertex_pair in self.body_vertices:
            v0 = vertex_pair[0] 
            v1 = vertex_pair[1]
            v0_transformed = T_SB_estim.dot(np.array([[v0[0]], [v0[1]], [v0[2]], [1]]))
            v1_transformed = T_SB_estim.dot(np.array([[v1[0]], [v1[1]], [v1[2]], [1]]))
            #self.axes.plot3D(*zip(v0_transformed[0:2], v1_transformed[0:2]), color="b")
            self.axes.plot([v0_transformed[0], v1_transformed[0]], [v0_transformed[1], v1_transformed[1]], color="g")

        # draw estimate landmarks
        for i in range(self.num_landmarks):       
            self.axes.scatter([self.m_estim[0][i]], [self.m_estim[1][i]], color=colormap(landmarks_quality[i]))  
            #self.axes.plot([self.ksi_estim[0][0], self.m_estim[0][i]], [self.ksi_estim[1][0], self.m_estim[1][i]], color="m", linestyle="dashed", linewidth=1.0)
            self.axes.plot([self.ksi_estim[0][0], self.ksi_estim[0][0] + 100.0*np.cos(z_mi[0][i]+self.ksi_estim[2][0])], [self.ksi_estim[1][0], self.ksi_estim[1][0] + 100.0*np.sin(z_mi[0][i]+self.ksi_estim[2][0])], color=colormap(landmarks_quality[i]), linestyle="dashed", linewidth=1.0)

        W, V = np.linalg.eig(self.S_estim[0 : 1 + 1 , 0 : 1 +1])
        idx = W.argsort()[::-1]   
        W = W[idx]
        V = V[:,idx]

        if W[0] < 0 or W[1] < 0:
            print("NEGATIVE EIGENVALUE (position)")

        # augment eigenvectors array to pseudo-3D
        V = np.concatenate((np.concatenate((V, np.array([[0.0, 0.0]])), axis=0), np.array([[0.0], [0.0], [0.0]])), axis=1)

        ellipsoid_t = np.linspace(0, 2*np.pi, 25)
        ellipsoid_x = np.sqrt(5.991*W[0]) * np.cos(ellipsoid_t)
        ellipsoid_y = np.sqrt(5.991*W[1]) * np.sin(ellipsoid_t)
        ellipsoid_x_transformed = ellipsoid_x.copy()
        ellipsoid_y_transformed = ellipsoid_y.copy()
        for j in range(len(ellipsoid_x_transformed)):
            T_estim_trans = np.concatenate((np.concatenate((axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(0.0)), np.array([[self.ksi_estim[0][0]], [self.ksi_estim[1][0]], [0]])), axis=1),
                                            np.array([[0, 0, 0, 1]])), axis=0)
            ellipsoid_i_xy = T_estim_trans.dot(np.concatenate((np.concatenate((V,np.array([[0],[0],[0]])), axis=1), np.array([[0,0,0,1]])), axis=0)).dot(np.array([[ellipsoid_x_transformed[j]], [ellipsoid_y_transformed[j]], [0], [1]]))
            ellipsoid_x_transformed[j] = ellipsoid_i_xy[0]
            ellipsoid_y_transformed[j] = ellipsoid_i_xy[1]
        self.axes.plot(ellipsoid_x_transformed, ellipsoid_y_transformed, color="c")

        for i in range(self.num_landmarks):
            W, V = np.linalg.eig(self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 3 + 2*i : 3 + 2*i+1 + 1])
            idx = W.argsort()[::-1]   
            W = W[idx]
            V = V[:,idx]      
      
            if W[0] < 0 or W[1] < 0:
                print("NEGATIVE EIGENVALUE (", i, ")")

            # augment eigenvectors array to pseudo-3D
            V = np.concatenate((np.concatenate((V, np.array([[0.0, 0.0]])), axis=0), np.array([[0.0], [0.0], [0.0]])), axis=1)

            ellipsoid_t = np.linspace(0, 2*np.pi, 25)
            ellipsoid_x = np.sqrt(5.991*W[0]) * np.cos(ellipsoid_t)
            ellipsoid_y = np.sqrt(5.991*W[1]) * np.sin(ellipsoid_t)
            ellipsoid_x_transformed = ellipsoid_x.copy()
            ellipsoid_y_transformed = ellipsoid_y.copy()
            for j in range(len(ellipsoid_x_transformed)):
                T_mi_estim_trans = np.concatenate((np.concatenate((axangles.axangle2mat(np.array([0, 0, 1]), np.deg2rad(0.0)), np.array([[self.m_estim[0][i]], [self.m_estim[1][i]], [0]])), axis=1),
                                                   np.array([[0, 0, 0, 1]])), axis=0)
                ellipsoid_i_xy = T_mi_estim_trans.dot(np.concatenate((np.concatenate((V,np.array([[0],[0],[0]])), axis=1), np.array([[0,0,0,1]])), axis=0)).dot(np.array([[ellipsoid_x_transformed[j]], [ellipsoid_y_transformed[j]], [0], [1]]))
                ellipsoid_x_transformed[j] = ellipsoid_i_xy[0]
                ellipsoid_y_transformed[j] = ellipsoid_i_xy[1]
            self.axes.plot(ellipsoid_x_transformed, ellipsoid_y_transformed, color=colormap(landmarks_quality[i]))
     
        # replace deprecated landmarks (odometry-like SLAM)
        print(*self.landmarks_quality_history)
        for i in range(self.num_landmarks):
            if self.landmarks_quality_history[i] >= self.landmarks_quality_counter_threshold:
                self.landmarks_quality_history[i] = -1

                self.m_groundtruth[:,[i]] = np.copy(self.ksi_groundtruth[0:1 + 1]) + np.array([[randomSign()*self.rng.uniform(1.0, 5.0)],[randomSign()*self.rng.uniform(1.0, 5.0)]])
                z_mi_dxy_init = self.m_groundtruth[:,[i]] - self.ksi_groundtruth[0:1 + 1]
                z_mi_dr_init = np.linalg.norm(z_mi_dxy_init) + self.rng.normal(0.0, 1.0)  #initialize at a depth of +-1.0m sigma (assumes we have some disparity-based delayed initialization pipeline for depth)
                z_mi_dtheta_init = np.arctan2(z_mi_dxy_init[1][0], z_mi_dxy_init[0][0]) - self.ksi_groundtruth[2][0] + np.max(np.array([-2.0*self.sigma_bearing, np.min(np.array([2.0*self.sigma_bearing, self.rng.normal(0.0, self.sigma_bearing)]))]))   #initialize based on sigma_bearing
                z_mi_dtheta_init = wrapToPi(z_mi_dtheta_init)
                z_mi[:,[i]] = np.array([[z_mi_dtheta_init]]) 

                self.m_estim[:,[i]] = np.array([[self.ksi_estim[0][0] + z_mi_dr_init*np.cos(self.ksi_estim[2][0] + z_mi_dtheta_init)], \
                                                [self.ksi_estim[1][0] + z_mi_dr_init*np.sin(self.ksi_estim[2][0] + z_mi_dtheta_init)]])

                self.S_estim[:                       , 3 + 2*i : 3 + 2*i+1 + 1] = np.zeros((3+2*self.num_landmarks, 2))
                self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , :]                       = np.zeros((2, 3+2*self.num_landmarks))
                self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 3 + 2*i : 3 + 2*i+1 + 1] = self.S_m_new

        # cache required values for next-iteration use (e.g. as previous hat robot pose estimate)
        self.ksi_hat_previous  = np.copy(ksi_hat)

        self.canvas.draw()

    # def on_update_values(self): 
    #     for i in range(self.m_groundtruth_text): 
    #         m_i_groundtruth_text = self.m_groundtruth_text[i].text().split(',')
    #         self.m_groundtruth[0][i] = float(m_i_groundtruth_text[0])
    #         self.m_groundtruth[1][i] = float(m_i_groundtruth_text[1])
    #     #print(self.m_groundtruth)

        #self.on_draw()

    def create_main_frame(self):
        self.main_frame = QWidget()

        self.dpi = 100
        self.fig = Figure((5.0, 15.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        self.axes = self.fig.add_subplot(111) #, projection='3d', proj_type='ortho'
        self.axes.set_aspect('equal')
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        

        self.m_groundtruth_text = list()
        for i in range(self.num_landmarks):
            self.m_groundtruth_text.append(QLineEdit())
            self.m_groundtruth_text[i].setMinimumWidth(65)
            self.m_groundtruth_text[i].setFixedWidth(65)
        #     self.connect(self.m_groundtruth_text[i], SIGNAL('editingFinished()'), self.on_update_values)

        hbox_rb = QHBoxLayout()
        for w in [ QLabel('landmarks x,y'), QLabel('#1'), self.m_groundtruth_text[0], QLabel('#2'), self.m_groundtruth_text[1], QLabel('#3'), self.m_groundtruth_text[2], QLabel('#4'), self.m_groundtruth_text[3], QLabel('#5'), self.m_groundtruth_text[4]]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_rb)
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Non-Holonomic Robot")
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

