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

from itertools import product, combinations

from numpy.random import default_rng

class Point2D:  
    def __init__(self, x, y):
        self.x = x
        self.y = y

    # Returns point -to- point a distance
    def distanceFromPoint(self, a):
        return np.sqrt(((a.x-self.x)**2)+((a.y-self.y)**2))

    # Returns point -to- line a-b distance
    def distanceFromLine(self, a, b):
        #np.abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))/np.sqrt(((x2-x1)**2)+((y2-y1)**2))
        return np.abs((b.x-a.x)*(a.y-self.y)-(a.x-self.x)*(b.y-a.y))/np.sqrt(((b.x-a.x)**2)+((b.y-a.y)**2))
    
    # Checks if point is between a-b
    def isBetween(self, a, b):
        crossproduct = (self.y - a.y) * (b.x - a.x) - (self.x - a.x) * (b.y - a.y)

        # compare versus epsilon for floating point values, or != 0 if using integers
        if abs(crossproduct) > 1e-9:
            return False

        dotproduct = (self.x - a.x) * (b.x - a.x) + (self.y - a.y)*(b.y - a.y)
        if dotproduct < 0:
            return False

        squaredlengthba = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y)
        if dotproduct > squaredlengthba:
            return False

        return True

    @staticmethod
    def linesIntersection(a, b, c, d):
        #intersect_x = ((x2*y1-x1*y2)*(x4-x3) - (x4*y3-x3*y4)*(x2-x1))/((x2-x1)*(y4-y3)-(x4-x3)*(y2-y1))
        #intersect_y = ((x2*y1-x1*y2)*(y4-y3) - (x4*y3-x3*y4)*(y2-y1))/((x2-x1)*(y4-y3)-(x4-x3)*(y2-y1))
        intersect_x = ((b.x*a.y-a.x*b.y)*(d.x-c.x) - (d.x*c.y-c.x*d.y)*(b.x-a.x))/((b.x-a.x)*(d.y-c.y)-(d.x-c.x)*(b.y-a.y))
        intersect_y = ((b.x*a.y-a.x*b.y)*(d.y-c.y) - (d.x*c.y-c.x*d.y)*(b.y-a.y))/((b.x-a.x)*(d.y-c.y)-(d.x-c.x)*(b.y-a.y))
        return Point2D(intersect_x, intersect_y)


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

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Non-Holonomic Robot')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        qApp.installEventFilter(self)

        self.fire_event = False

        self.rng = default_rng()

        self.m_groundtruth = np.transpose(np.array([[3.0, 2.0], \
                                                   [5.0, -1.0], \
                                                   [-5.0, -1.0], \
                                                   [-4.0, -5.0], \
                                                   [-3.0, 5.0]]))
        self.m_0_groundtruth.setText(str(self.m_groundtruth[0][0])+","+str(self.m_groundtruth[1][0]))
        self.m_1_groundtruth.setText(str(self.m_groundtruth[0][1])+","+str(self.m_groundtruth[1][1]))
        self.m_2_groundtruth.setText(str(self.m_groundtruth[0][2])+","+str(self.m_groundtruth[1][2]))
        self.m_3_groundtruth.setText(str(self.m_groundtruth[0][3])+","+str(self.m_groundtruth[1][3]))
        self.m_4_groundtruth.setText(str(self.m_groundtruth[0][4])+","+str(self.m_groundtruth[1][4]))

        self.num_landmarks = len(np.transpose(self.m_groundtruth))

        self.sigma_Rwheel = 1.0 * np.deg2rad(5.0)
        self.sigma_Lwheel = 1.0 * np.deg2rad(5.0)

        self.sigma_range = 2.0 * 0.25
        self.sigma_angle = 2.0 * np.deg2rad(2.5)
   
        self.Q = 1.0 * np.array([[0.1,               0], \
                                 [0,   np.deg2rad(1.0)]])
        self.R = 1.0 * np.array([[0.25,               0], \
                                 [0,    np.deg2rad(2.5)]])

        self.m_estim = np.copy(self.m_groundtruth)
        # for i in range(0, len(self.m_groundtruth)):
        #     self.m_estim[i] = 5 * [self.m_groundtruth[i][0], self.m_groundtruth[i][1]] / np.linalg.norm(self.m_groundtruth[i])  #initialize at 5m range

        ksi_0 = np.array([[0.0], \
                          [0.0], \
                          [np.deg2rad(0.0)]])  # [x, y, theta]
        self.ksi_groundtruth = ksi_0
        self.ksi_estim = ksi_0

        self.S_ksi_estim = 1.0 * np.array([[0.001,     0,     0], \
                                           [0,     0.001,     0], \
                                           [0,         0, 0.001]])
        self.S_m_new = 1.0 * np.array([[0.5,   0], \
                                       [0,   0.1]])
        self.S_estim = np.concatenate((np.concatenate((self.S_ksi_estim, np.zeros((3, self.m_groundtruth.size))), axis=1), \
                                       np.concatenate((np.transpose(np.zeros((3, self.m_groundtruth.size))), np.kron(np.eye(self.num_landmarks,dtype=float),self.S_m_new)), axis=1)), axis=0)

        self.u_t = np.array([[0.0], \
                             [0.0]])  # input forward and turn velocity

        self.z_mi = np.empty((2,self.num_landmarks))  # measurements vector for all landmarks
        self.z_mi[:,:] = np.NaN

        self.timer_period = 0.10
        self.timer_value = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_timer)
        self.timer.start(self.timer_period * 1000)

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

        #self.on_draw()

    def eventFilter(self, source, event):
        if event.type() == QEvent.KeyPress:
            #print('KeyPress: %s [%r]' % (event.key(), source))
            if event.key() == Qt.Key_Up:  # up
                self.u_t[0][0] = self.u_t[0][0] + 0.5
            elif event.key() == 16777237:  # down
                self.u_t[0][0] = self.u_t[0][0] - 0.5
            elif event.key() == 16777234:  # left
                self.u_t[1][0] = self.u_t[1][0] + 0.1
            elif event.key() == 16777236:  # right
                self.u_t[1][0] = self.u_t[1][0] - 0.1
            elif event.key() == Qt.Key_Space:  # spacebar
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

        if np.abs(self.u_t[0][0]) > 0:
            if np.abs(self.u_t[0][0]) <= 0.05:
                self.u_t[0][0] = 0
            else:
                self.u_t[0][0] = self.u_t[0][0] - np.sign(self.u_t[0][0]) * 0.05
        if np.abs(self.u_t[0][0]) > 0:
            if np.abs(self.u_t[1][0]) <= 0.01:
                self.u_t[1][0] = 0
            else:
                self.u_t[1][0] = self.u_t[1][0] - np.sign(self.u_t[1][0]) * 0.01

        # update low-level Joint-Space commands (wheels) & add noise
        u_R = self.u_t[0][0] + self.u_t[1][0] * self.half_width + np.max(np.array([-2.0*self.sigma_Rwheel, np.min(np.array([2.0*self.sigma_Rwheel, self.rng.standard_normal()*self.sigma_Rwheel]))]))
        u_L = self.u_t[0][0] - self.u_t[1][0] * self.half_width + np.max(np.array([-2.0*self.sigma_Lwheel, np.min(np.array([2.0*self.sigma_Lwheel, self.rng.standard_normal()*self.sigma_Lwheel]))]))

        # calculate Body Frame kinematics based on joint-space velocities
        u_vel = (u_R + u_L)/2;
        u_rot = (u_R - u_L)/(2 * self.half_width)

        x_dot = u_vel * np.cos( self.ksi_groundtruth[2][0] )
        y_dot = u_vel * np.sin( self.ksi_groundtruth[2][0] )
        theta_dot = u_rot;  

        # update Space Frame kinematics 
        self.ksi_groundtruth[0][0] = self.ksi_groundtruth[0][0] + x_dot * self.timer_period
        self.ksi_groundtruth[1][0] = self.ksi_groundtruth[1][0] + y_dot * self.timer_period
        self.ksi_groundtruth[2][0] = self.ksi_groundtruth[2][0] + theta_dot * self.timer_period

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

        scale_full = 1.0
        self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), color=['r'])
        self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), color=['g'])
        self.axes.set_xlim(-10.0, 10.0)
        self.axes.set_ylim(-10.0, 10.0)

        # draw groundtruth body
        for vertex_pair in self.body_vertices:
            v0 = vertex_pair[0] 
            v1 = vertex_pair[1]
            v0_transformed = T_SB.dot(np.array([[v0[0]], [v0[1]], [v0[2]], [1]]))
            v1_transformed = T_SB.dot(np.array([[v1[0]], [v1[1]], [v1[2]], [1]]))
            #self.axes.plot3D(*zip(v0_transformed[0:2], v1_transformed[0:2]), color="b")
            self.axes.plot([v0_transformed[0], v1_transformed[0]], [v0_transformed[1], v1_transformed[1]], color="b")

        # draw groundtruth landmarks
        for i in range(self.num_landmarks):
            self.axes.scatter([self.m_groundtruth[0][i]], [self.m_groundtruth[1][i]], color="k")
            self.axes.plot([self.ksi_groundtruth[0][0], self.m_groundtruth[0][i]], [self.ksi_groundtruth[1][0], self.m_groundtruth[1][i]], color="k", linestyle="dashed", linewidth=1.0)

        # form range bearing measurements, corrupted by noise
        z_mi = np.empty((2,self.num_landmarks)) 
        z_mi[:,:] = np.NaN
        for i in range(self.num_landmarks):  
            z_mi_dx_groundtruth = self.m_groundtruth[0][i] - self.ksi_groundtruth[0][0]
            z_mi_dy_groundtruth = self.m_groundtruth[1][i] - self.ksi_groundtruth[1][0]
            z_mi_dr_groundtruth = np.linalg.norm(np.array([z_mi_dx_groundtruth,z_mi_dy_groundtruth]))
            z_mi[0][i] = z_mi_dr_groundtruth + np.max(np.array([-2.0*self.sigma_range, np.min(np.array([2.0*self.sigma_range, self.rng.standard_normal()*self.sigma_range]))]))
            z_mi[1][i] = np.arctan2(z_mi_dy_groundtruth, z_mi_dx_groundtruth) - self.ksi_groundtruth[2][0] + np.max(np.array([-2.0*self.sigma_angle, np.min(np.array([2.0*self.sigma_angle, self.rng.standard_normal()*self.sigma_angle]))]))

        # predicted estimate
        ksi_hat = np.copy(self.ksi_estim)
        ksi_hat[0][0] = ksi_hat[0][0] + self.timer_period * self.u_t[0][0]*np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)  #(-self.u_t[0][0]/self.u_t[1][0])*np.sin(self.ksi_estim[2][0]) + (self.u_t[0][0]/self.u_t[1][0])*np.sin(self.ksi_estim[2][0] + self.u_t[1][0]*self.timer_period);
        ksi_hat[1][0] = ksi_hat[1][0] + self.timer_period * self.u_t[0][0]*np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)  #(self.u_t[0][0]/self.u_t[1][0])*np.cos(self.ksi_estim[2][0]) + (-self.u_t[0][0]/self.u_t[1][0])*np.cos(self.ksi_estim[2][0] + self.u_t[1][0]*self.timer_period);
        ksi_hat[2][0] = ksi_hat[2][0] + self.timer_period * self.u_t[1][0]

        F_ksi = np.array([[1, 0, self.timer_period * -self.u_t[0][0]*np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                          [0, 1, self.timer_period * self.u_t[0][0]*np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                          [0, 0, 1]])
        F_u = np.array([[self.timer_period * np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period), self.timer_period * self.u_t[0][0]*-0.5*self.timer_period*np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                        [self.timer_period * np.sin(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period), self.timer_period * self.u_t[0][0]*0.5*self.timer_period*np.cos(self.ksi_estim[2][0] + 0.5*self.u_t[1][0]*self.timer_period)], \
                        [0                                                                                      , self.timer_period]])
        S_ksi_hat = F_ksi.dot(self.S_estim[0:3,0:3]).dot(np.transpose(F_ksi)) + F_u.dot(self.Q).dot(np.transpose(F_u))

        # form range bearing predictions and Jacobian entries
        z_mi_hat = np.empty((2,self.num_landmarks)) 
        z_mi_hat[:,:] = np.NaN
        H_mi = np.empty((2 * self.num_landmarks, 5)) 
        H_mi[:,:] = np.NaN
        for i in range(self.num_landmarks):  
            z_mi_dx_hat = self.m_estim[0][i] - self.ksi_groundtruth[0][0]
            z_mi_dy_hat = self.m_estim[1][i] - self.ksi_groundtruth[1][0]
            z_mi_dr_hat = np.linalg.norm(np.array([z_mi_dx_hat,z_mi_dy_hat]))
            z_mi_hat[0][i] = z_mi_dr_hat
            z_mi_hat[1][i] = np.arctan2(z_mi_dy_hat, z_mi_dx_hat) - self.ksi_groundtruth[2][0]
            
            H_mi[2*i:2*i+1 + 1,:] = np.array([[-z_mi_dx_hat/z_mi_dr_hat    , -z_mi_dy_hat/z_mi_dr_hat     , 0 ,  z_mi_dx_hat/z_mi_dr_hat     , z_mi_dy_hat/z_mi_dr_hat], \
                                              [z_mi_dy_hat/(z_mi_dr_hat**2), -z_mi_dx_hat/(z_mi_dr_hat**2), -1, -z_mi_dy_hat/(z_mi_dr_hat**2), z_mi_dx_hat/(z_mi_dr_hat**2)]]) 

        # Landmark i update
        for i in range(self.num_landmarks):
            S_hat = np.zeros((5,5))
            S_hat[0 : 2 + 1   , 0 : 2 + 1]   = np.copy(S_ksi_hat)
            S_hat[3 : 3+1 + 1 , 3 : 3+1 + 1] = np.copy(self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 3 + 2*i : 3 + 2*i+1 + 1])
            S_hat[0 : 2 + 1   , 3 : 3+1 + 1] = np.copy(self.S_estim[0               : 2 + 1 , 3 + 2*i : 3+ 2*i+1 + 1])
            S_hat[3 : 3+1 + 1 , 0 : 2 + 1]   = np.copy(self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 0       : 2 + 1])
            K_mi = S_hat.dot(np.transpose(H_mi[2*i : 2*i+1 + 1 , :])).dot(np.linalg.inv(H_mi[2*i : 2*i+1 + 1 , :].dot(S_hat).dot(np.transpose(H_mi[2*i : 2*i+1 + 1 , :])) + self.R))

            x_hat = np.concatenate((ksi_hat,
                                    self.m_estim[:,[i]]), axis=0)
            x_upd = x_hat + K_mi.dot( z_mi[:,[i]] - z_mi_hat[:,[i]] )

            S_upd = ( np.eye(3+2) - K_mi.dot(H_mi[2*i : 2*i+1 + 1 , :]) ).dot(S_hat);

            self.ksi_estim    = np.copy(x_upd[0   : 2 + 1,[0]])
            self.m_estim[:,[i]] = np.copy(x_upd[2+1 : 2+2 + 1][:])
            self.S_estim[0       : 2 + 1         , 0       : 2 + 1]         = np.copy(S_upd[0 : 2 + 1   , 0 : 2 + 1])
            self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 3 + 2*i : 3 + 2*i+1 + 1] = np.copy(S_upd[3 : 3+1 + 1 , 3 : 3+1 + 1])
            self.S_estim[0       : 2 + 1         , 3 + 2*i : 3 + 2*i+1 + 1] = np.copy(S_upd[0 : 2 + 1   , 3 : 3+1 + 1])
            self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 0       : 2 + 1]         = np.copy(S_upd[3 : 3+1 + 1 , 0 : 2 + 1])

        # Update pose covariance (auxilliary) variable
        self.S_ksi_estim = np.copy(self.S_estim[0:2 + 1,0:2 + 1])

        # draw estimate
        T_SB_estim = np.concatenate((np.concatenate((axangles.axangle2mat(np.array([0, 0, 1]), self.ksi_estim[2][0]), np.array([[self.ksi_estim[0][0]], [self.ksi_estim[1][0]], [0]])), axis=1),
                                     np.array([[0, 0, 0, 1]])), axis=0)

        self.quiver_Bp_estim = T_SB_estim.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
        self.quiver_Bx_estim = T_SB_estim.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
        self.quiver_By_estim = T_SB_estim.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
        self.quiver_Bz_estim = T_SB_estim.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

        scale_full = 1.0
        self.axes.quiver(self.quiver_Bp_estim[0], self.quiver_Bp_estim[1], scale_full*(self.quiver_Bx_estim[0]-self.quiver_Bp_estim[0]), scale_full*(self.quiver_Bx_estim[1]-self.quiver_Bp_estim[1]), color=['r'])
        self.axes.quiver(self.quiver_Bp_estim[0], self.quiver_Bp_estim[1], scale_full*(self.quiver_By_estim[0]-self.quiver_Bp_estim[0]), scale_full*(self.quiver_By_estim[1]-self.quiver_Bp_estim[1]), color=['g'])
  
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
            self.axes.scatter([self.m_estim[0][i]], [self.m_estim[1][i]], color="m")  
            self.axes.plot([self.ksi_estim[0][0], self.m_estim[0][i]], [self.ksi_estim[1][0], self.m_estim[1][i]], color="m", linestyle="dashed", linewidth=1.0)

        W, V = np.linalg.eig(self.S_ksi_estim)
        idx = W.argsort()[::-1]   
        W = W[idx]
        V = V[:,idx]

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
        self.axes.plot(ellipsoid_x_transformed, ellipsoid_y_transformed, color="m")

        for i in range(self.num_landmarks):
            W, V = np.linalg.eig(self.S_estim[3 + 2*i : 3 + 2*i+1 + 1 , 3 + 2*i : 3 + 2*i+1 + 1])
            idx = W.argsort()[::-1]   
            W = W[idx]
            V = V[:,idx]      
      
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
            self.axes.plot(ellipsoid_x_transformed, ellipsoid_y_transformed, color="m")
     
        self.canvas.draw()

    def on_update_values(self):       
        m_0_groundtruth = self.m_0_groundtruth.text().split(',')
        self.m_groundtruth[0][0] = float(m_0_groundtruth[0])
        self.m_groundtruth[1][0] = float(m_0_groundtruth[1])
        m_1_groundtruth = self.m_1_groundtruth.text().split(',')
        self.m_groundtruth[0][1] = float(m_1_groundtruth[0])
        self.m_groundtruth[1][1] = float(m_1_groundtruth[1])
        m_2_groundtruth = self.m_2_groundtruth.text().split(',')
        self.m_groundtruth[0][2] = float(m_2_groundtruth[0])
        self.m_groundtruth[1][2] = float(m_2_groundtruth[1])
        m_3_groundtruth = self.m_3_groundtruth.text().split(',')
        self.m_groundtruth[0][3] = float(m_3_groundtruth[0])
        self.m_groundtruth[1][3] = float(m_3_groundtruth[1])
        m_4_groundtruth = self.m_4_groundtruth.text().split(',')
        self.m_groundtruth[0][4] = float(m_4_groundtruth[0])
        self.m_groundtruth[1][4] = float(m_4_groundtruth[1])
        #print(self.m_groundtruth)

        self.val_vel = self.sld_vel.value()
        self.val_lookahead = self.sld_lookahead.value()

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
        

        self.m_0_groundtruth = QLineEdit()
        self.m_0_groundtruth.setMinimumWidth(65)
        self.m_0_groundtruth.setFixedWidth(65)
        self.connect(self.m_0_groundtruth, SIGNAL('editingFinished()'), self.on_update_values)

        self.m_1_groundtruth = QLineEdit()
        self.m_1_groundtruth.setMinimumWidth(65)
        self.m_1_groundtruth.setFixedWidth(65)
        self.connect(self.m_1_groundtruth, SIGNAL('editingFinished()'), self.on_update_values)

        self.m_2_groundtruth = QLineEdit()
        self.m_2_groundtruth.setMinimumWidth(65)
        self.m_2_groundtruth.setFixedWidth(65)
        self.connect(self.m_2_groundtruth, SIGNAL('editingFinished()'), self.on_update_values)

        self.m_3_groundtruth = QLineEdit()
        self.m_3_groundtruth.setMinimumWidth(65)
        self.m_3_groundtruth.setFixedWidth(65)
        self.connect(self.m_3_groundtruth, SIGNAL('editingFinished()'), self.on_update_values)

        self.m_4_groundtruth = QLineEdit()
        self.m_4_groundtruth.setMinimumWidth(65)
        self.m_4_groundtruth.setFixedWidth(65)
        self.connect(self.m_4_groundtruth, SIGNAL('editingFinished()'), self.on_update_values)

        self.sld_vel = DoubleSlider(Qt.Horizontal)
        self.sld_vel.setMinimum(0.0)
        self.sld_vel.setMaximum(2.5)
        self.sld_vel.setValue(0.0)
        self.sld_vel.setTracking(True)
        self.sld_vel.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_vel, SIGNAL('valueChanged(int)'), self.on_update_values)

        self.sld_lookahead = DoubleSlider(Qt.Horizontal)
        self.sld_lookahead.setMinimum(0.25)
        self.sld_lookahead.setMaximum(10.0)
        self.sld_lookahead.setValue(5.0)
        self.sld_lookahead.setTracking(True)
        self.sld_lookahead.setTickPosition(QSlider.TicksBelow)
        self.connect(self.sld_lookahead, SIGNAL('valueChanged(int)'), self.on_update_values)

        # self.rb_g0 = QCheckBox('Guidance #1')
        # self.rb_g0.setChecked(True)
        # self.connect(self.rb_g0, SIGNAL('stateChanged(int)'), self.on_update_values)

        hbox_vel = QHBoxLayout()
        for w in [ QLabel('vel'), QLabel('0'), self.sld_vel, QLabel('2.5')]:
            hbox_vel.addWidget(w)
            hbox_vel.setAlignment(w, Qt.AlignVCenter)

        hbox_lookahead = QHBoxLayout()
        for w in [ QLabel('lookahead'), QLabel('0.25'), self.sld_lookahead, QLabel('10')]:
            hbox_lookahead.addWidget(w)
            hbox_lookahead.setAlignment(w, Qt.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ QLabel('landmarks x,y'), QLabel('#1'), self.m_0_groundtruth, QLabel('#2'), self.m_1_groundtruth, QLabel('#3'), self.m_2_groundtruth, QLabel('#4'), self.m_3_groundtruth, QLabel('#5'), self.m_4_groundtruth]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_vel)
        vbox.addLayout(hbox_lookahead)
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

