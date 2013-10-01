#coding=utf8
# virtualenv test_py2
import sys
import time
import pandas
import numpy as np
from colors import ColorWheel
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import logging
from scipy.integrate import ode, odeint

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARN)
handler = logging.StreamHandler()
formatter = logging.Formatter('[%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

qt_app = QtGui.QApplication(sys.argv)

class WorldSimple(object):
    def __init__(self,birthrate=0.03,deathrate=0.01,regenerationrate=0.1,
                 burdenrate=0.02,economyaim=1,growthrate=0.05,resolution=1e5):
        self.birthrate = birthrate
        self.deathrate = deathrate
        self.regenerationrate = regenerationrate
        self.burdenrate = burdenrate
        self.economyaim = economyaim
        self.growthrate = growthrate
        self.resolution = resolution
        self.time = np.linspace(0,250,resolution)

    def dx(self,x,time):
        population,burden,economy = x
        quality = burden**(-1)
        birth = self.birthrate * population * quality * economy
        death = population * self.deathrate * burden
        ecocide = self.burdenrate * economy * population
        if quality > 1:
            regeneration = self.regenerationrate * burden
        else:
            regeneration = self.regenerationrate
        economicgrowth = self.growthrate * economy * burden \
                * (1-(economy*burden)/self.economyaim)
        return np.array([birth-death,ecocide-regeneration,economicgrowth])

    def _solve(self,x0=[1,1,1]):
        return odeint(self.dx,x0,self.time,
                      full_output=True,printmessg=False,mxhnil=0)

    def solve(self):
        res,info = self._solve()
        if info['message'] == "Integration successful.":
            res = pandas.DataFrame(res,
                columns=['population','burden','economy'])
            res['time'] = self.time
            return res

class DataFramePlot(object):
    def __init__(self,plotwidget):
        self.plots = {}
        self.colorwheel = ColorWheel()
        self.plotwidget = plotwidget
        self.legend = plotwidget.addLegend()


    def addItem(self,dataframe,column):
        self.plots[column] = self.plotwidget.plotItem.plot(
            dataframe['time'],dataframe[column].values)
        #,symbol='-',
            #pen=tuple(self.colorwheel.next().rgb),name=column
            #)
        #self.plotwidget.addItem(self.plots[column])

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
                    self.plots[column].setData(dataframe.index,dataframe[column])

class Parameter(QtGui.QGroupBox):
    def __init__(self,name,value=1,xmin=0,xmax=10,step=1.0):
        super(Parameter,self).__init__(name)
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
        self.textbox = QtGui.QLineEdit()
        self.textbox.setMaxLength(30)
        self.textbox.setFixedWidth(60)
        self.textbox.setText(str(value))
        self.textbox.setValidator(QtGui.QIntValidator(0,neededsteps))
        layout.addWidget(self.textbox)
        layout.addWidget(self.slider)
        self.setLayout(layout)
        self.slider.valueChanged.connect(self.on_value_slider_changed)
        self.textbox.textChanged.connect(self.on_value_textbox_changed)

    valueChanged = QtCore.pyqtSignal(str,float)

    def on_value_slider_changed(self,value):
        self.value = self.step*value+self.min
        self.textbox.setText(str(self.value))
        self.valueChanged.emit(self.param_name,self.value)

    def on_value_textbox_changed(self,value):
        if value:
            self.value = float(value)
            self.slider.setValue((self.value-self.min)/self.step)
            self.valueChanged.emit(self.param_name,self.value)

class WorldSimpleGui(object):
    """ segementation fault bug, see:
        https://groups.google.com/d/msg/pyqtgraph/juoqAiXABYY/BYpvCQeWSgMJ
    """
    def __init__(self):
        self.plotwin = pg.PlotWindow(title="Simple Limits of Growth World Model")
        self.plotwin.resize(1000,600)
        self.world = WorldSimple(resolution=1000)
        # initial solution and plotting
        self.dataframeplot = DataFramePlot(self.plotwin)
        self.dataframeplot.create_plots(self.world.solve())
        self.controllerwin = QtGui.QWidget()
        self.controllerwin.setWindowTitle('Parameters')
        layout = QtGui.QVBoxLayout()
        for param in ('birthrate','deathrate','regenerationrate',
                      'burdenrate','economyaim','growthrate'):
            slider = Parameter(param,step=0.01)
            slider.valueChanged.connect(self.change)
            layout.addWidget(slider)
        self.controllerwin.setLayout(layout)

    def run(self):
        self.plotwin.show()
        self.controllerwin.show()
        qt_app.exec_()

    def change(self,param, value):
        oldvalue = getattr(self.world,str(param))
        setattr(self.world,str(param),value)
        t1 = time.time()
        res = self.world.solve()
        if res:
            self.dataframeplot.update_plots(res)
            logger.info('recalculated in {0}s'.format(time.time()-t1))
        else:
            logger.warn('failed with parameter "{0}" equal {1}'.format(param,value))

if __name__ == '__main__':
    m = WorldSimpleGui()
    m.run()
