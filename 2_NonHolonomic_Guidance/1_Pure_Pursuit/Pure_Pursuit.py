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


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Non-Holonomic Robot')

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()


        self.val_vel = 0
        self.val_lookahead = 5

        self.trajectory = list() 
        self.trajectory.append([5, 5])
        self.trajectory.append([7.5, -7.5])
        self.trajectory.append([-7.5, -5])
        self.trajectory.append([-5, 5])
        self.trajectory.append([0, 7.5])
        self.traj_xy_1.setText(str(self.trajectory[0][0])+","+str(self.trajectory[0][1]))
        self.traj_xy_2.setText(str(self.trajectory[1][0])+","+str(self.trajectory[1][1]))
        self.traj_xy_3.setText(str(self.trajectory[2][0])+","+str(self.trajectory[2][1]))
        self.traj_xy_4.setText(str(self.trajectory[3][0])+","+str(self.trajectory[3][1]))
        self.traj_xy_5.setText(str(self.trajectory[4][0])+","+str(self.trajectory[4][1]))

        self.traj_active_segment = 0 
        self.lookahead_intersect_point_previous = [np.NAN, np.NAN]
        self.lookahead_intersect_point = [np.NAN, np.NAN]

        self.ksi_groundtruth = [0, 0, np.deg2rad(0)]  # [x, y, theta]

        self.timer_period = 0.10
        self.timer_value = 0
        self.timer = QTimer()
        self.timer.timeout.connect(self.on_timer)
        self.timer.start(self.timer_period * 1000)

        self.scenario = 0

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
        half_length = 0.25
        half_width = 0.15
        half_height = 0.05
        self.body_vertices = list()
        self.body_vertices.append([[-half_length, -half_width, -half_height], [-half_length, -half_width, half_height]])
        self.body_vertices.append([[-half_length, -half_width, -half_height], [half_length, -half_width, -half_height]])
        self.body_vertices.append([[-half_length, -half_width, -half_height], [half_length, -half_width, -half_height]])
        self.body_vertices.append([[-half_length, -half_width, half_height], [-half_length, half_width, half_height]])
        self.body_vertices.append([[-half_length, -half_width, half_height], [half_length, -half_width, half_height]])
        self.body_vertices.append([[-half_length, half_width, -half_height], [-half_length, half_width, half_height]])
        self.body_vertices.append([[-half_length, half_width, -half_height], [half_length, half_width, -half_height]])
        self.body_vertices.append([[-half_length, half_width, half_height], [half_length, half_width, half_height]])
        self.body_vertices.append([[half_length, -half_width, -half_height], [half_length, -half_width, half_height]])
        self.body_vertices.append([[half_length, -half_width, -half_height], [half_length, half_width, -half_height]])
        self.body_vertices.append([[half_length, -half_width, half_height], [half_length, half_width, half_height]])
        self.body_vertices.append([[half_length, half_width, -half_height], [half_length, half_width, half_height]])

        self.half_width = half_width

        # point-line distance to choose closest segment to initialize active_segment
        min_dist = 10000.0
        for i in range(len(self.trajectory)):
            a = Point2D(self.ksi_groundtruth[0], self.ksi_groundtruth[1])
            b = Point2D(self.trajectory[i][0], self.trajectory[i][1])
            if i < len(self.trajectory)-1:
                c = Point2D(self.trajectory[i+1][0], self.trajectory[i+1][1])
            else:
                c = Point2D(self.trajectory[0][0], self.trajectory[0][1])
            dist = a.distanceFromLine(b, c)
            if dist < min_dist:
                min_dist = dist
                self.traj_active_segment = i

        #self.on_draw()

    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def steer(self):
        # iterate different angles around the vehicle, find minimum-distance intersection, and look-ahead distance intersection (if it exists)
        min_lookahead_difference = 10000.0
        self.lookahead_intersect_point_previous = self.lookahead_intersect_point
        for d_a in range(0, 360):
            a = Point2D(self.ksi_groundtruth[0], self.ksi_groundtruth[1])
            b = Point2D(a.x + 1*np.cos(np.deg2rad(d_a)), a.y + 1*np.sin(np.deg2rad(d_a)))

            c = Point2D(self.trajectory[self.traj_active_segment][0], self.trajectory[self.traj_active_segment][1])
            if self.traj_active_segment < len(self.trajectory)-1:
                d = Point2D(self.trajectory[self.traj_active_segment+1][0], self.trajectory[self.traj_active_segment+1][1])
            else:
                d = Point2D(self.trajectory[0][0], self.trajectory[0][1])

            intersect = Point2D.linesIntersection(a, b, c, d)
            position = Point2D(self.ksi_groundtruth[0], self.ksi_groundtruth[1])
            intersect_length = intersect.distanceFromPoint(position)

            lookahead_difference = np.abs(self.val_lookahead-intersect_length)
            if lookahead_difference < min_lookahead_difference and intersect_length <= self.val_lookahead:

                if self.traj_active_segment < len(self.trajectory)-1:
                    next_waypoint = Point2D(self.trajectory[self.traj_active_segment+1][0], self.trajectory[self.traj_active_segment+1][1])
                else:
                    next_waypoint = Point2D(self.trajectory[0][0], self.trajectory[0][1])

                if (next_waypoint.distanceFromPoint(intersect) < \
                    next_waypoint.distanceFromPoint(Point2D(self.lookahead_intersect_point_previous[0], self.lookahead_intersect_point_previous[1])) ) or \
                   (np.isnan(self.lookahead_intersect_point_previous[0]) or np.isnan(self.lookahead_intersect_point_previous[1])):
                    # drop ambiguous solutions by having to move closer to the next waypoint
                    min_lookahead_difference = lookahead_difference
                    self.lookahead_intersect_point = [intersect.x, intersect.y]

        if not np.isnan(self.lookahead_intersect_point[0]) and not np.isnan(self.lookahead_intersect_point[1]):
            # valid lookahead intersection point, steer towards that
            intersect_point = self.lookahead_intersect_point
            omega_command = (2 * self.val_vel * np.sin(np.arctan2(intersect_point[1]-self.ksi_groundtruth[1], intersect_point[0]-self.ksi_groundtruth[0])-self.ksi_groundtruth[2])) / self.val_lookahead

            # check if outside trajectory line segment endpoints to advance self.traj_active_segment
            position = Point2D(self.ksi_groundtruth[0], self.ksi_groundtruth[1])
            intersection = Point2D(self.lookahead_intersect_point[0], self.lookahead_intersect_point[1])
            a = Point2D(self.trajectory[self.traj_active_segment][0], self.trajectory[self.traj_active_segment][1])
            if self.traj_active_segment < len(self.trajectory)-1:
                b = Point2D(self.trajectory[self.traj_active_segment+1][0], self.trajectory[self.traj_active_segment+1][1])
            else:
                b = Point2D(self.trajectory[0][0], self.trajectory[0][1])
            if (not intersection.isBetween(a, b) and position.distanceFromPoint(b) < self.val_lookahead) or \
               (position.distanceFromPoint(b) < self.val_lookahead/2):
                self.traj_active_segment = self.traj_active_segment + 1
                if self.traj_active_segment >= len(self.trajectory):
                    self.traj_active_segment = 0
            #print(intersection.isBetween(b, a))
        else:
            # not valid, just steer towards active segment's final waypoint
            if self.traj_active_segment < len(self.trajectory)-1:
                intersect_point = [self.trajectory[self.traj_active_segment+1][0], self.trajectory[self.traj_active_segment+1][1]]
            else:
                intersect_point = [self.trajectory[0][0], self.trajectory[0][1]]
            omega_command = (2 * self.val_vel * np.sin(np.arctan2(intersect_point[1]-self.ksi_groundtruth[1], intersect_point[0]-self.ksi_groundtruth[0])-self.ksi_groundtruth[2])) / self.val_lookahead

        return omega_command   

    def on_timer(self): 
        self.timer_value = self.timer_value + self.timer_period

        if np.abs(self.val_vel) < 1e-3:
            val_omega = 0
        else:
            val_omega = self.steer()

        # update low-level Joint-Space commands (wheels)
        u_right = self.val_vel + val_omega * self.half_width
        u_left = self.val_vel - val_omega * self.half_width
   
        # calculate Body Frame kinematics based on joint-space velocities
        vel = (u_right + u_left)/2;
        omega = (u_right - u_left)/(2 * self.half_width);
        
        x_dot = vel * np.cos( self.ksi_groundtruth[2] );
        y_dot = vel * np.sin( self.ksi_groundtruth[2] );
        theta_dot = omega;

        # update Space Frame kinematics 
        self.ksi_groundtruth[0] = self.ksi_groundtruth[0] + x_dot * self.timer_period;
        self.ksi_groundtruth[1] = self.ksi_groundtruth[1] + y_dot * self.timer_period;
        self.ksi_groundtruth[2] = self.ksi_groundtruth[2] + theta_dot * self.timer_period;

        #print(self.timer_value)

        self.on_draw()

    def on_draw(self):
        self.axes.clear()        
        self.axes.grid(True)

        if self.scenario == 0:

            T_SB = np.concatenate((np.concatenate((axangles.axangle2mat(np.array([0, 0, 1]), self.ksi_groundtruth[2]), np.array([[self.ksi_groundtruth[0]], [self.ksi_groundtruth[1]], [0]])), axis=1),
                                   np.array([[0, 0, 0, 1]])), axis=0)

            self.quiver_Bp = T_SB.dot(np.concatenate((self.S_p, np.array([[1]])), axis=0))
            self.quiver_Bx = T_SB.dot(np.concatenate((self.quiver_Sx, np.array([[1]])), axis=0))
            self.quiver_By = T_SB.dot(np.concatenate((self.quiver_Sy, np.array([[1]])), axis=0))
            self.quiver_Bz = T_SB.dot(np.concatenate((self.quiver_Sz, np.array([[1]])), axis=0))

            # these are just to scale arrows of different coordinate systems to better distinguish between them
            scale_small = 0.25
            scale_full = 1.0
            self.axes.quiver(self.S_p[0], self.S_p[1], scale_small*self.quiver_Sx[0], scale_small*self.quiver_Sx[1], color=['r'])
            self.axes.quiver(self.S_p[0], self.S_p[1], scale_small*self.quiver_Sy[0], scale_small*self.quiver_Sy[1], color=['g'])
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], scale_full*(self.quiver_Bx[0]-self.quiver_Bp[0]), scale_full*(self.quiver_Bx[1]-self.quiver_Bp[1]), color=['r'])
            self.axes.quiver(self.quiver_Bp[0], self.quiver_Bp[1], scale_full*(self.quiver_By[0]-self.quiver_Bp[0]), scale_full*(self.quiver_By[1]-self.quiver_Bp[1]), color=['g'])
            self.axes.set_xlim(-10.0, 10.0)
            self.axes.set_ylim(-10.0, 10.0)

            for vertex_pair in self.body_vertices:
                v0 = vertex_pair[0] 
                v1 = vertex_pair[1]
                v0_transformed = T_SB.dot(np.array([[v0[0]], [v0[1]], [v0[2]], [1]]))
                v1_transformed = T_SB.dot(np.array([[v1[0]], [v1[1]], [v1[2]], [1]]))
                #self.axes.plot3D(*zip(v0_transformed[0:2], v1_transformed[0:2]), color="b")
                self.axes.plot([v0_transformed[0], v1_transformed[0]], [v0_transformed[1], v1_transformed[1]], color="b")

            for i in range(len(self.trajectory)-1):
                self.axes.plot([self.trajectory[i][0], self.trajectory[i+1][0]], [self.trajectory[i][1], self.trajectory[i+1][1]], color="k")
            self.axes.plot([self.trajectory[len(self.trajectory)-1][0], self.trajectory[0][0]], [self.trajectory[len(self.trajectory)-1][1], self.trajectory[0][1]], color="k")

            if not np.isnan(self.lookahead_intersect_point[0]) and not np.isnan(self.lookahead_intersect_point[1]):
                self.axes.plot([self.ksi_groundtruth[0], self.lookahead_intersect_point[0]], [self.ksi_groundtruth[1], self.lookahead_intersect_point[1]], color="r")

        self.canvas.draw()

    def on_update_values(self):       
        traj_xy_1 = self.traj_xy_1.text().split(',')
        self.trajectory[0][0] = float(traj_xy_1[0])
        self.trajectory[0][1] = float(traj_xy_1[1])
        traj_xy_2 = self.traj_xy_2.text().split(',')
        self.trajectory[1][0] = float(traj_xy_2[0])
        self.trajectory[1][1] = float(traj_xy_2[1])
        traj_xy_3 = self.traj_xy_3.text().split(',')
        self.trajectory[2][0] = float(traj_xy_3[0])
        self.trajectory[2][1] = float(traj_xy_3[1])
        traj_xy_4 = self.traj_xy_4.text().split(',')
        self.trajectory[3][0] = float(traj_xy_4[0])
        self.trajectory[3][1] = float(traj_xy_4[1])
        traj_xy_5 = self.traj_xy_5.text().split(',')
        self.trajectory[4][0] = float(traj_xy_5[0])
        self.trajectory[4][1] = float(traj_xy_5[1])
        #print(self.trajectory)

        self.val_vel = self.sld_vel.value()
        self.val_lookahead = self.sld_lookahead.value()

        #if self.rb_t0.isChecked() and self.scenario != 0:
        #    self.scenario = 0

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
        

        self.traj_xy_1 = QLineEdit()
        self.traj_xy_1.setMinimumWidth(50)
        self.traj_xy_1.setFixedWidth(50)
        self.connect(self.traj_xy_1, SIGNAL('editingFinished()'), self.on_update_values)

        self.traj_xy_2 = QLineEdit()
        self.traj_xy_2.setMinimumWidth(50)
        self.traj_xy_2.setFixedWidth(50)
        self.connect(self.traj_xy_2, SIGNAL('editingFinished()'), self.on_update_values)

        self.traj_xy_3 = QLineEdit()
        self.traj_xy_3.setMinimumWidth(50)
        self.traj_xy_3.setFixedWidth(50)
        self.connect(self.traj_xy_3, SIGNAL('editingFinished()'), self.on_update_values)

        self.traj_xy_4 = QLineEdit()
        self.traj_xy_4.setMinimumWidth(50)
        self.traj_xy_4.setFixedWidth(50)
        self.connect(self.traj_xy_4, SIGNAL('editingFinished()'), self.on_update_values)

        self.traj_xy_5 = QLineEdit()
        self.traj_xy_5.setMinimumWidth(50)
        self.traj_xy_5.setFixedWidth(50)
        self.connect(self.traj_xy_5, SIGNAL('editingFinished()'), self.on_update_values)

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

        self.rb_g0 = QCheckBox('Guidance #1')
        self.rb_g0.setChecked(True)
        self.connect(self.rb_g0, SIGNAL('stateChanged(int)'), self.on_update_values)

        hbox_vel = QHBoxLayout()
        for w in [ QLabel('vel'), QLabel('0'), self.sld_vel, QLabel('2.5')]:
            hbox_vel.addWidget(w)
            hbox_vel.setAlignment(w, Qt.AlignVCenter)

        hbox_lookahead = QHBoxLayout()
        for w in [ QLabel('lookahead'), QLabel('0.25'), self.sld_lookahead, QLabel('10')]:
            hbox_lookahead.addWidget(w)
            hbox_lookahead.setAlignment(w, Qt.AlignVCenter)

        hbox_rb = QHBoxLayout()
        for w in [ self.rb_g0, QLabel('trajectory x,y'), QLabel('#1'), self.traj_xy_1, QLabel('#2'), self.traj_xy_2, QLabel('#3'), self.traj_xy_3, QLabel('#4'), self.traj_xy_4, QLabel('#5'), self.traj_xy_5]:
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

