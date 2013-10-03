#coding=utf8
# virtualenv test_py2
import os
import sys
import time
import pandas
import numpy as np
from colors import ColorWheel
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import logging
from scipy.integrate import ode, odeint
import gettext
import locale
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode
from pkg_resources import resource_filename

logger = logging.getLogger(__name__)

file_formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
file_handler = logging.FileHandler('limitsofgrowth.log')
file_handler.setFormatter(file_formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

console_formatter = logging.Formatter('[%(levelname)s] %(message)s')
console_handler = logging.StreamHandler()
console_handler.setFormatter(console_formatter)
console_handler.setLevel(logging.WARN)
logger.addHandler(console_handler)

try:
    t = gettext.translation('limitsofgrowth',resource_filename(__name__,'data/locale/'))
except:
    t = gettext.NullTranslations()

_ = t.ugettext


qt_app = QtGui.QApplication(sys.argv)


class WorldSimple(object):
    params = {'birthrate':{'initial':0.03,'max':1,'min':0},
              'deathrate':{'initial':0.01,'max':1,'min':0},
              'regenerationrate':{'initial':0.1,'max':1,'min':0},
              'burdenrate':{'initial':0.02,'max':1,'min':0},
              'economyaim':{'initial':10,'max':100,'min':0},
              'growthrate':{'initial':0.05,'max':1,'min':0}}

    cols = ['population','burden','economy']
    def __init__(self,resolution=1e3):
        self.time = np.linspace(0,250,resolution)
        self.sensitivities = []
        for column in self.cols:
            for param in self.params:
                self.sensitivities.append('{0},{1}'.format(column,param))

    def create_dx(self,params):
        birthrate = params['birthrate']
        deathrate = params['deathrate']
        regenerationrate = params['regenerationrate']
        burdenrate = params['burdenrate']
        economyaim = params['economyaim']
        growthrate = params['growthrate']
        def func(time,x):
            """
            the change of the system variables population, burden and economy
            x: [population,burden,economy]
            """
            quality = x[1]**(-1)
            birth = birthrate * x[0] * quality * x[2]
            death = x[0] * deathrate * x[1]
            ecocide = burdenrate * x[2] * x[0]
            if quality > 1:
                regeneration = regenerationrate * x[1]
            else:
                regeneration = regenerationrate
            economicgrowth = growthrate * x[2] * x[1] \
                    * (1-(x[2]*x[1])/economyaim)
            return [birth-death,ecocide-regeneration,economicgrowth]
        return func
    @staticmethod
    def dx_with_parameters(time,x,p):
        """
        the change of the system variables population, burden and economy
        x: [population,burden,economy]
        """
        quality = x[1]**(-1)
        birth = p[0] * x[0] * quality * x[2]
        death = x[0] * p[1] * x[1]
        ecocide = p[3] * x[2] * x[0]
        if quality > 1:
            regeneration = p[2] * x[1]
        else:
            regeneration = p[2]
        economicgrowth = p[5] * x[2] * x[1] \
                * (1-(x[2]*x[1])/p[4])
        return [birth-death,ecocide-regeneration,economicgrowth]

    def x_odeint(self,params):
        func = self.create_dx(params)
        res,info = odeint(lambda x,t:func(t,x),[1.0,1.0,1.0], self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe = pandas.DataFrame(res,
                columns=['population','burden','economy'])
            dataframe['time'] = self.time
            return dataframe

    def x_cvode(self,params):
        problem = Explicit_Problem(self.create_dx(params), [1.0,1.0,1.0],0)
        sim = CVode(problem)
        t,x = sim.simulate(250,1000)
        dataframe = pandas.DataFrame(x,
                columns=['population','burden','economy'])
        dataframe['time'] = t
        return dataframe

    def create_ds(self,params,x):
        birthrate = params['birthrate']
        deathrate = params['deathrate']
        regenerationrate = params['regenerationrate']
        burdenrate = params['burdenrate']
        economyaim = params['economyaim']
        growthrate = params['growthrate']
        def func(time,s1d):
            """
            the change of sensitivities of the system variables
            """
            # since the ode solver only work with 1 dimensional arrays, lets just
            # internaly work with multidimensional arrays
            # get nearest time index
            index = (x['time'] - time).abs().argmin()
            population = x['population'][index]
            burden = x['burden'][index]
            economy = x['economy'][index]
            s = s1d.reshape((3,6))
            fnachx = np.array([
                [birthrate/burden*economy-deathrate*burden,
                birthrate*population*(-1)/burden**2*economy-population*deathrate,
                birthrate*population/burden],
                [burdenrate*economy,
                regenerationrate if burden < 1 else 0,
                burdenrate*population],
                [0,
                growthrate*economy-2*growthrate*economy**2*burden/economyaim,
                growthrate*burden-2*growthrate*burden**2*economy/economyaim]
            ])

            fnachp = np.array([
                [population*burden*economy,
                - population+burden,0,0,0,0],
                [0,0,-burden if burden<1 else -1,
                economy*population,0,0],
                [0,0,0,0,growthrate*economy**2*burden**2/economyaim**2,
                economyaim*burden*(1-(economy*burden/economyaim))]
            ])
            return (fnachx.dot(s)+fnachp).reshape((18,))
        return func

    def s_odeint(self,params):
        x = self.x_odeint(params)
        s0 = np.zeros(18)
        func = self.create_ds(params,x)
        res,info = odeint(lambda s,t:func(t,s),np.zeros(18), self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe = pandas.DataFrame(res,columns=self.sensitivities)
            dataframe['time'] = self.time
            return dataframe

    def s_cvode(self,params):
        x = self.x_odeint(params)
        s0 = np.zeros(18)
        problem = Explicit_Problem(self.create_ds(params,x),s0,0)
        sim = CVode(problem)
        t,s = sim.simulate(250,1000)
        dataframe = pandas.DataFrame(s,columns=self.sensitivities)
        dataframe['time'] = t
        return dataframe

    def s_cvode_natural(self,params):
        problem = Explicit_Problem(self.dx_with_parameters, [1.0,1.0,1.0],0,(
            params['birthrate'],params['deathrate'],params['regenerationrate'],
            params['burdenrate'],params['economyaim'],params['growthrate']))
        sim = CVode(problem)
        sim.report_continuously = True
        t,x = sim.simulate(250,self.time.shape[0]-1)
        dataframe = pandas.DataFrame(x,
                columns=['population','burden','economy'])
        dataframe['time'] = self.time
        d = {}
        sens = np.array(sim.p_sol)
        for i,col in enumerate(self.cols):
            for j,param in enumerate(
                ('birthrate','deathrate','regenerationrate',
                 'burdenrate','economyaim','growthrate')):
                d['{0},{1}'.format(col,param)] = sens[j,:,i]
        dataframe_sens = pandas.DataFrame(d,index=self.time)
        dataframe_sens['time'] = self.time
        return dataframe,dataframe_sens

    def s_ode(self,params):
        x0 = np.ones(3)
        s0 = np.zeros(18)
        solver_x = ode(self.dx).set_integrator('dopri5')
        solver_x.set_initial_value(x0,0).set_f_params(params,)

        solver_s = ode(self.ds).set_integrator('dopri5')
        solver_s.set_initial_value(s0,0).set_f_params(params,x0)

        dt = 1
        t1 = 250
        sensitivities = []
        for column in self.cols:
            for param in self.params:
                sensitivities.append('{0},{1}'.format(column,param))
        sol = pandas.DataFrame(np.empty((t1/dt,22)),columns=self.cols+sensitivities+['time'])
        i = 0
        #return solver_x,solver_s
        while solver_x.successful() and solver_s.successful() and solver_x.t < t1:
            solver_x.integrate(solver_x.t+dt)
            sol.iloc[i][self.cols] = solver_x.y
            sol.iloc[i]['time'] = solver_x.t
            solver_s.set_f_params(params,solver_x.y)
            solver_s.integrate(solver_s.t+dt)
            sol.iloc[i][sensitivities] = solver_s.y
            i += 1
        return sol

class DataFramePlot(object):
    def __init__(self,plotwidget):
        self.plots = {}
        self.colorwheel = ColorWheel(start=0.6)
        self.plotwidget = plotwidget
        self.legend = plotwidget.addLegend()

    def addItem(self,dataframe,column):
        self.plots[column] = self.plotwidget.plotItem.plot(
            dataframe.index.values,dataframe[column].values,
            pen=tuple(self.colorwheel.next().rgb),name=_(column)
        )

    def create_plots(self,dataframe):
        for column in dataframe:
            if column != 'time':
                self.addItem(dataframe,column)

    def update_plots(self,dataframe):
        for column in dataframe:
            if column != 'time':
                if not column in self.plots:
                    self.addItem(dataframe,column)
                else:
                    self.plots[column].setData(dataframe['time'],dataframe[column])

class Parameter(QtGui.QGroupBox):
    def __init__(self,name,value=1,xmin=0,xmax=10,step=1.0):
        super(Parameter,self).__init__(_(name))
        self.param_name = name
        self.value = value
        layout = QtGui.QHBoxLayout()
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal, parent=self)
        neededsteps = (xmax-xmin)/step
        self.max = xmax
        self.min = xmin
        self.step = step
        self.slider.setRange(0,neededsteps)
        self.slider.setValue((value-xmin)/step)  # set the initial position
        self.slider.setPageStep(neededsteps/10)
        self.slider.setTickInterval(neededsteps/10)
        self.slider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.textbox = QtGui.QLineEdit()
        self.textbox.setMaxLength(30)
        self.textbox.setFixedWidth(60)
        self.textbox.setText(str(value))
        self.textbox.setValidator(QtGui.QDoubleValidator(0.0,float(neededsteps),3))
        layout.addWidget(self.textbox)
        layout.addWidget(self.slider)
        self.setLayout(layout)
        self.slider.valueChanged.connect(self.on_value_slider_changed)
        self.textbox.editingFinished.connect(self.on_value_textbox_changed)

    valueChanged = QtCore.pyqtSignal(str,float)

    def on_value_slider_changed(self,value):
        self.value = self.step*value+self.min
        self.textbox.setText(str(self.value))
        self.valueChanged.emit(self.param_name,self.value)

    def on_value_textbox_changed(self):
        value = str(self.textbox.text())
        if value:
            self.value = float(value)
            self.slider.setValue((self.value-self.min)/self.step)
            self.valueChanged.emit(self.param_name,self.value)

class MatplotlibPlot(QtGui.QWidget):
    def __init__(self):
        super(MatplotlibPlot,self).__init__()
        self.resize(1000,600)
        self.create_main_frame()

    def on_draw(self):
        """ Redraws the figure
        """
        # clear the axes and redraw the plot anew
        #
        self.axes.clear()
        for column in self.dataframe:
            self.axes.plot(self.dataframe.index,self.dataframe[column].values,label=column)
        self.axes.legend()
        self.canvas.draw()

    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points

        QtGui.QMessageBox.information(self, "Click!", msg)

    def create_main_frame(self):

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)

        # Bind the 'pick' event for clicking on one of the bars
        #
        self.canvas.mpl_connect('pick_event', self.on_pick)

        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        self.setLayout(vbox)

    def add_additional_stuff(self):
        # Other GUI controls
        #
        self.textbox = QtGui.QLineEdit()
        self.textbox.setMinimumWidth(200)
        self.textbox.editingFinished.connect(self.on_draw)

        self.draw_button = QtGui.QPushButton("&Draw")
        self.draw_button.clicked.connect(self.on_draw)

        self.grid_cb = QtGui.QCheckBox("Show &Grid")
        self.grid_cb.setChecked(False)
        self.grid_cb.stateChanged.connect(self.on_draw)

        slider_label = QtGui.QLabel('Bar width (%):')
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slider.setRange(1, 100)
        self.slider.setValue(20)
        self.slider.setTracking(True)
        self.slider.setTickPosition(QtGui.QSlider.TicksBothSides)
        self.slider.valueChanged.connect(self.on_draw)

        #
        # Layout with box sizers
        #
        hbox = QtGui.QHBoxLayout()

        for w in [  self.textbox, self.draw_button, self.grid_cb,
                    slider_label, self.slider]:
            hbox.addWidget(w)
            hbox.setAlignment(w, QtCore.Qt.AlignVCenter)

        vbox = self.layout()
        vbox.addLayout(hbox)


class WorldSimpleGui(pg.PlotWidget):
    """ segementation fault bug, see:
        https://groups.google.com/d/msg/pyqtgraph/juoqAiXABYY/BYpvCQeWSgMJ
    """
    def __init__(self):
        super(WorldSimpleGui,self).__init__()
        self.world = WorldSimple(resolution=1000)
        self.create_main_window()
        self.create_fullscreen_shortcut()
        self.create_parameter_widget()
        self.create_main_controls_widget()
        self.create_initial_plot()

    def create_initial_plot(self):
        self.parameters = {param:value['initial']
                            for param,value in self.world.params.items()}
        self.dataframeplot = DataFramePlot(self)
        self.dataframeplot.create_plots(self.world.x_odeint(self.parameters))

    def create_main_window(self):
        self.win = QtGui.QMainWindow()
        self.win.setCentralWidget(self)
        self.win.setWindowTitle("Simple Limits of Growth World Model")
        self.win.resize(1000,600)
        self.tools = QtGui.QTabWidget()
        dock = QtGui.QDockWidget('Controllers')
        dock.setWidget(self.tools)
        self.win.addDockWidget(QtCore.Qt.RightDockWidgetArea,dock)

    def create_fullscreen_shortcut(self):
        self.shcut_fullscreen = QtGui.QShortcut(self.win)
        self.shcut_fullscreen.setKey("F11")
        self.shcut_fullscreen.activated.connect(self.toogleFullscreen)

    def create_main_controls_widget(self):
        widget = QtGui.QWidget()
        widget.resize(300,widget.height())
        layout = QtGui.QVBoxLayout()
        sens_button = QtGui.QPushButton('Sensitivity Analysis')
        sens_button.clicked.connect(self.plot_sensitivity)
        layout.addWidget(sens_button)
        widget.setLayout(layout)
        self.tools.addTab(widget,'Tasks')

    def create_parameter_widget(self):
        widget = QtGui.QWidget()
        widget.resize(300,widget.height())
        layout = QtGui.QVBoxLayout()
        for param,value in self.world.params.items():
            slider = Parameter(param,value=value['initial'],
                               xmin=value['min'],xmax=value['max'],step=0.01)
            slider.valueChanged.connect(self.change)
            layout.addWidget(slider)
        widget.setLayout(layout)
        self.tools.addTab(widget,'Parameters')

    def plot_sensitivity(self):
        x,s = self.world.s_cvode_natural(self.parameters)
        data = s[['{0},economyaim'.format(x) for x in self.world.cols]]
        self.p = MatplotlibPlot()
        self.p.dataframe = data
        self.p.on_draw()
        self.p.show()

    def toogleFullscreen(self):
        if self.win.isFullScreen():
            self.win.showNormal()
        else:
            self.win.showFullScreen()

    def run(self):
        self.win.show()
        #self.parameter_widget.show()
        qt_app.exec_()

    def change(self,param, value):
        oldvalue = self.parameters[str(param)]
        self.parameters[str(param)] = value
        t1 = time.time()
        res = self.world.x_odeint(self.parameters)
        if res:
            self.dataframeplot.update_plots(res)
            logger.info('recalculated in {0}s'.format(time.time()-t1))
        else:
            logger.warn('failed with parameter "{0}" equal {1}'.format(param,value))


def main():
    m = WorldSimpleGui()
    m.run()

if __name__ == '__main__':
    main()
