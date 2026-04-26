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
from scipy.spatial import KDTree

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


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


def wrapToPi(radians):
    while radians >= np.pi:
        radians = radians - 2.0*np.pi
    while radians < -np.pi:
        radians = radians + 2.0*np.pi
    return radians

def setupBody(x, y, width, height, theta):
    half_width = width/2
    half_height = height/2

    body_vertices = list()
    body_vertices.append([-half_width, -half_height])
    body_vertices.append([-half_width, +half_height])
    body_vertices.append([+half_width, +half_height])
    body_vertices.append([+half_width, -half_height])

    for i, vertex in enumerate(body_vertices):
        rotated_vertex = np.array([[np.cos(theta), -np.sin(theta)], \
                                   [np.sin(theta),  np.cos(theta)]]).dot(vertex)
        body_vertices[i] = np.array([x, y]) + rotated_vertex

    body_edges = list()
    body_edges.append([[body_vertices[0][0], body_vertices[0][1]], [body_vertices[1][0], body_vertices[1][1]]])
    body_edges.append([[body_vertices[1][0], body_vertices[1][1]], [body_vertices[2][0], body_vertices[2][1]]])
    body_edges.append([[body_vertices[2][0], body_vertices[2][1]], [body_vertices[3][0], body_vertices[3][1]]])
    body_edges.append([[body_vertices[3][0], body_vertices[3][1]], [body_vertices[0][0], body_vertices[0][1]]])

    body_polygon = Polygon([(body_vertices[0][0], body_vertices[0][1]), (body_vertices[1][0], body_vertices[1][1]), (body_vertices[2][0], body_vertices[2][1]), (body_vertices[3][0], body_vertices[3][1])])

    return body_edges, body_polygon

def steer_to_goal(pose_curr, pose_goal, discretization):
    trajectory_to_goal = []

    # exact steering: no kinodynamic constraints
    theta_at_goal = np.arctan2(pose_goal[1]-pose_curr[1],pose_goal[0]-pose_curr[0])
    pose_at_goal = np.array([pose_goal[0], pose_goal[1], theta_at_goal])

    pose_delta_xy = pose_at_goal[0:2] - pose_curr[0:2]
    pose_delta_xy_length = np.linalg.norm(pose_delta_xy)

    pose_xy_step_length = discretization
    while pose_xy_step_length < pose_delta_xy_length:
        pose_step_xy = pose_curr[0:2] + pose_xy_step_length*(pose_delta_xy/pose_delta_xy_length)
        pose_step_theta = pose_curr[2]
        trajectory_to_goal.append(np.array([pose_step_xy[0], pose_step_xy[1], pose_step_theta]))
        pose_xy_step_length = pose_xy_step_length + discretization
    
    trajectory_to_goal.append(pose_at_goal)

    trajectory_to_goal_numposes = len(trajectory_to_goal)
    for i, trajectory_pose in enumerate(trajectory_to_goal):
        delta_theta_total = trajectory_to_goal[-1][2] - trajectory_to_goal[0][2]
        delta_theta_i = i*(delta_theta_total/trajectory_to_goal_numposes)
        trajectory_to_goal[i][2] = trajectory_to_goal[0][2] + delta_theta_i

    return trajectory_to_goal

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

class Node:
    def __init__(self, value, parent):
        self.value = value
        self.parent = parent

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Non-Holonomic Robot')

        self.obstacles = np.transpose(np.array([[5.0, 2.0, 5.0, 7.0, np.deg2rad(45.0)], \
                                                [3.0, -6.0, 4.0, 3.0, np.deg2rad(30.0)], \
                                                [-6.0, 1.0, 7.5, 2.0, np.deg2rad(15.0)], \
                                                [-4.0, -5.0, 2.0, 6.0, np.deg2rad(-45.0)], \
                                                [-3.0, 6.5, 3.0, 5.0, np.deg2rad(-25.0)]]))
        self.num_obstacles = len(np.transpose(self.obstacles))

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        qApp.installEventFilter(self)

        self.fire_event = False

        self.rng = default_rng()

        for i in range(self.num_obstacles): 
            self.obstacles_text[i].setText(str(self.obstacles[0][i])+","+str(self.obstacles[1][i]))

        self.ksi_init = np.array([0.0, \
                                  0.0, \
                                  np.deg2rad(0.0)])  # [x, y, theta]

        self.ksi_sample = np.array([0.0, \
                                    0.0])  # [x, y]
        
        self.ksi_steer = np.copy(self.ksi_init)

        self.root = Node(self.ksi_init[0:2], None)

        self.ksi_goal = np.array([7.5, \
                                  -7.5])  # [x, y]

        self.sample_goal_counter = 0

        self.goal_reached_thres = 1.0

        self.max_conn_length = 3.5

        self.kdtree_data = [self.root]

        self.timer_period = 0.01
        self.timer = QTimer()
        self.timer.timeout.connect(self.plan)
        self.timer.start(self.timer_period * 1000)
        self.timer_timeout_counter = 0

        self.half_width = 10.0
        self.half_height = 10.0

        self.robot_width = 1.5
        self.robot_height = 0.75
        self.robot_edges_init, self.robot_polygon_init = setupBody(0.0, 0.0, self.robot_width, self.robot_height, np.deg2rad(0.0))

        robot_edges, robot_polygon = setupBody(0.0, 0.0, self.robot_width, self.robot_height, np.deg2rad(0.0))
        self.trajectory_robot_edges = [robot_edges]
        self.trajectory_robot_polygons = [robot_polygon]

        self.obstacle_edges = []
        self.obstacle_polygons = []
        for i in range(self.num_obstacles): 
            obstacle_edges, obstacle_polygon = setupBody(self.obstacles[0][i], self.obstacles[1][i], self.obstacles[2][i], self.obstacles[3][i], self.obstacles[4][i])
            self.obstacle_edges.append(obstacle_edges)
            self.obstacle_polygons.append(obstacle_polygon)

        self.tree_edges_vis = []

        self.tree_path_to_goal_vis = []

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

    def eventFilter(self, source, event):
        if event.type() == QEvent.KeyPress:
            #print('KeyPress: %s [%r]' % (event.key(), source))
            if event.key() == Qt.Key_Space:  # spacebar
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
    
    def plan(self):

        self.sample_goal_counter = self.sample_goal_counter + 1
        if self.sample_goal_counter >= 10:
            # sample the goal itself every so often
            self.sample_goal_counter = 0
            self.ksi_sample = np.copy(self.ksi_goal)
        else:
            # randomly sample and constrain within workspace dimensions
            self.ksi_sample[:] = np.array([self.rng.normal(0.0, self.half_width),self.rng.normal(0.0, self.half_height)])
            self.ksi_sample[0] = np.clip(self.ksi_sample[0], -self.half_width+1.1*self.robot_width/2, self.half_width-1.1*self.robot_width/2)
            self.ksi_sample[1] = np.clip(self.ksi_sample[1], -self.half_height+1.1*self.robot_height/2, self.half_height-1.1*self.robot_height/2)

        # find 1st nearest neighbor
        kdtree = KDTree([node.value[0:2] for node in self.kdtree_data])
        nn_distance, nn_index = kdtree.query(self.ksi_sample[0:2], k=1)
        nn_node = self.kdtree_data[nn_index]
      
        # adjust to limit single-connection step within reasonable max radius
        if nn_distance > self.max_conn_length:
            self.ksi_sample[0] = nn_node.value[0] + self.max_conn_length * (self.ksi_sample[0]-nn_node.value[0])/nn_distance
            self.ksi_sample[1] = nn_node.value[1] + self.max_conn_length * (self.ksi_sample[1]-nn_node.value[1])/nn_distance

        # steer to sample
        self.ksi_steer_trajectory = steer_to_goal(np.array([nn_node.value[0], nn_node.value[1], 0.0]), \
                                                  np.array([self.ksi_sample[0], self.ksi_sample[1], 0.0]), \
                                                  1.0)
        self.ksi_steer_goal = np.copy(self.ksi_steer_trajectory[-1])

        self.trajectory_robot_edges = []
        self.trajectory_robot_polygons = []
        for ksi_steer in self.ksi_steer_trajectory:
            robot_edges, robot_polygon = setupBody(ksi_steer[0], ksi_steer[1], self.robot_width, self.robot_height, ksi_steer[2])
            self.trajectory_robot_edges.append(robot_edges)
            self.trajectory_robot_polygons.append(robot_polygon)

        # check intersection with obstacles
        collision_free = True
        for obstacle_polygon in self.obstacle_polygons:
            for robot_polygon in self.trajectory_robot_polygons:
                if obstacle_polygon.intersects(robot_polygon):
                    collision_free = False
                    break

        if collision_free:
            newnode = Node(self.ksi_steer_goal, nn_node)
            self.kdtree_data.append(newnode)
            self.tree_edges_vis.append([nn_node.value, newnode.value])

            if np.linalg.norm(self.ksi_goal - self.ksi_steer_goal[0:2]) <= self.goal_reached_thres:
                currnode = newnode
                while not currnode.parent is None:
                    self.tree_path_to_goal_vis.append([currnode.value, currnode.parent.value])
                    currnode = currnode.parent

        self.timer_timeout_counter = self.timer_timeout_counter + 1
        if self.timer_timeout_counter >= 0:
            self.timer_timeout_counter = 0
            self.on_draw()
 
    def on_draw(self):
        self.axes.clear()        
        self.axes.grid(True)

        self.axes.set_xlim(-self.half_width, self.half_width)
        self.axes.set_ylim(-self.half_height, self.half_height)

        # draw initial config
        for edge in self.robot_edges_init:
            vertex_pair = edge
            v0 = vertex_pair[0]
            v1 = vertex_pair[1]
            self.axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color="k")

        # draw obstacles
        for i in range(self.num_obstacles):
            for edge in self.obstacle_edges[i]:
                vertex_pair = edge
                v0 = vertex_pair[0] 
                v1 = vertex_pair[1]
                self.axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color="b")

        # draw goal
        self.axes.scatter(self.ksi_goal[0], self.ksi_goal[1], color="r", s=500.0)  

        # draw sample config
        self.axes.scatter(self.ksi_sample[0], self.ksi_sample[1], color="m")  

        # draw steer-to-goal trajectory configurations
        for robot_edges in self.trajectory_robot_edges:
            for edge in robot_edges:
                vertex_pair = edge
                v0 = vertex_pair[0]
                v1 = vertex_pair[1]
                self.axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color="g")

        # draw tree edges
        for edge in self.tree_edges_vis:
            vertex_pair = edge
            v0 = vertex_pair[0]
            v1 = vertex_pair[1]
            self.axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color="y")

        # draw solution
        for edge in self.tree_path_to_goal_vis:
            vertex_pair = edge
            v0 = vertex_pair[0]
            v1 = vertex_pair[1]
            self.axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color="r")

        self.canvas.draw()

    def create_main_frame(self):
        self.main_frame = QWidget()

        self.dpi = 100
        self.fig = Figure((5.0, 15.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        self.axes = self.fig.add_subplot(111) #, projection='3d', proj_type='ortho'
        self.axes.set_aspect('equal')
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        

        self.obstacles_text = list()
        for i in range(self.num_obstacles):
            self.obstacles_text.append(QLineEdit())
            self.obstacles_text[i].setMinimumWidth(65)
            self.obstacles_text[i].setFixedWidth(65)
        #     self.connect(self.obstacles_text[i], SIGNAL('editingFinished()'), self.on_update_values)

        hbox_rb = QHBoxLayout()
        for w in [ QLabel('obstacles x,y'), QLabel('#1'), self.obstacles_text[0], QLabel('#2'), self.obstacles_text[1], QLabel('#3'), self.obstacles_text[2], QLabel('#4'), self.obstacles_text[3], QLabel('#5'), self.obstacles_text[4]]:
            hbox_rb.addWidget(w)
            hbox_rb.setAlignment(w, Qt.AlignVCenter)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_rb)
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Rapidly-exploring Random Trees")
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

