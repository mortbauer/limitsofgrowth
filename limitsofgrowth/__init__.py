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
    params = {'birthrate':{'initial':0.03,'max':1,'min':0},
              'deathrate':{'initial':0.01,'max':1,'min':0},
              'regenerationrate':{'initial':0.1,'max':1,'min':0},
              'burdenrate':{'initial':0.02,'max':1,'min':0},
              'economyaim':{'initial':10,'max':100,'min':0},
              'growthrate':{'initial':0.05,'max':1,'min':0}}

    def __init__(self,resolution=1e5):
        self.time = np.linspace(0,250,resolution)

    @staticmethod
    def dx(x,time,params):
        population,burden,economy = x
        quality = burden**(-1)
        birth = params['birthrate'] * population * quality * economy
        death = population * params['deathrate'] * burden
        ecocide = params['burdenrate'] * economy * population
        if quality > 1:
            regeneration = params['regenerationrate'] * burden
        else:
            regeneration = params['regenerationrate']
        economicgrowth = params['growthrate'] * economy * burden \
                * (1-(economy*burden)/params['economyaim'])
        return np.array([birth-death,ecocide-regeneration,economicgrowth])


    def solve(self,params):
        res,info = odeint(self.dx,[1.0,1.0,1.0],self.time,args=(params,),
                      full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            res = pandas.DataFrame(res,
                columns=['population','burden','economy'])
            res['time'] = self.time
            return res

class DataFramePlot(object):
    def __init__(self,plotwidget):
        self.plots = {}
        self.colorwheel = ColorWheel(start=0.6)
        self.plotwidget = plotwidget
        self.legend = plotwidget.addLegend()


    def addItem(self,dataframe,column):
        self.plots[column] = self.plotwidget.plotItem.plot(
            dataframe['time'].values,dataframe[column].values,
            pen=tuple(self.colorwheel.next().rgb),name=column
        )
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
                    self.plots[column].setData(dataframe['time'],dataframe[column])

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
        self.slider.setTickInterval(neededsteps/10)
        self.slider.setTickPosition(QtGui.QSlider.TicksBelow)
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

class WorldSimpleGui(pg.PlotWidget):
    """ segementation fault bug, see:
        https://groups.google.com/d/msg/pyqtgraph/juoqAiXABYY/BYpvCQeWSgMJ
    """
    def __init__(self):
        super(WorldSimpleGui,self).__init__()
        self.win = QtGui.QMainWindow()
        self.win.setCentralWidget(self)
        self.win.setWindowTitle("Simple Limits of Growth World Model")
        self.win.closeEvent = self.closeEvent
        self.win.resize(1000,600)
        self.world = WorldSimple(resolution=1000)
        # initial solution and plotting
        self.controllerwin = QtGui.QWidget()
        self.controllerwin.resize(300,self.controllerwin.height())
        self.controllerwin.move(1100,0)
        self.controllerwin.closeEvent = self.closeEvent
        self.controllerwin.setWindowTitle('Parameters')
        self.shcut1 = QtGui.QShortcut(self.win)
        self.shcut1.setKey("F11")
        self.connect(self.shcut1, QtCore.SIGNAL("activated()"), self.toogleFullscreen)

        self.params = {}
        layout = QtGui.QVBoxLayout()
        for param in self.world.params:
            d = self.world.params[param]
            self.params[param] = d['initial']
            slider = Parameter(param,value=d['initial'],
                               xmin=d['min'],xmax=d['max'],step=0.01)
            slider.valueChanged.connect(self.change)
            layout.addWidget(slider)
        self.dataframeplot = DataFramePlot(self)
        self.dataframeplot.create_plots(self.world.solve(self.params))
        self.controllerwin.setLayout(layout)

    def toogleFullscreen(self):
        if self.win.isFullScreen():
            self.win.showNormal()
        else:
            self.win.showFullScreen()

    def closeEvent(self,event):
        self.controllerwin.close()
        self.win.close()

    def run(self):
        self.win.show()
        self.controllerwin.show()
        qt_app.exec_()

    def change(self,param, value):
        oldvalue = self.params[str(param)]
        self.params[str(param)] = value
        t1 = time.time()
        res = self.world.solve(self.params)
        if res:
            self.dataframeplot.update_plots(res)
            logger.info('recalculated in {0}s'.format(time.time()-t1))
        else:
            logger.warn('failed with parameter "{0}" equal {1}'.format(param,value))

if __name__ == '__main__':
    m = WorldSimpleGui()
    m.run()
