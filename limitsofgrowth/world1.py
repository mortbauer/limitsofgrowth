#coding=utf8
# virtualenv test_py2
import pandas
import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import logging
from scipy.integrate import ode, odeint
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
formatter = logging.Formatter('[%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

# state space representation

# y: population, burden, economy
# x: birth, death, ecocide, regeneration,  economicgrowth
# u: birthrate, deathrate, regenerationrate, burdenrate, economyaim, growthrate

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
            self.regenerationrate
        economicgrowth = self.growthrate * economy * burden \
                * (1-(economy*burden)/self.economyaim)
        return np.array([birth-death,ecocide-regeneration,economicgrowth])

    def solve(self,x0=[1,1,1]):
        res = pandas.DataFrame(odeint(self.dx,x0,self.time),
            columns=['population','burden','economy'])
        res['time'] = self.time
        return res

class WorldSimpleGui(object):
    def __init__(self):
        self.world = WorldSimple()
        # set up the app and widgets
        win = pg.GraphicsWindow(title="very simple World Model")
        layout = QtGui.QGridLayout()
        win.setLayout(layout)
        text = QtGui.QLineEdit('enter text')
        layout.addWidget(text, 0, 1)
        plot = pg.PlotWidget()
        layout.addWidget(plot, 0, 0)
        # initial solution and plotting
        self._r = r = self.world.solve()
        self.population_plot = pg.PlotCurveItem(
            r['time'],r['population'],pen=(255,0,0),name='population')
        self.burden_plot = pg.PlotCurveItem(
            r['time'],r['burden'],pen=(0,255,0),name='burden')
        self.economy_plot = pg.PlotCurveItem(
            r['time'],r['economy'],pen=(0,0,255),name='economy')
        plot.addItem(self.population_plot)
        plot.addItem(self.burden_plot)
        plot.addItem(self.economy_plot)
        #self.layout.addWidget(self.plot, 0, 1, 3, 1)
        # set upt the parameter tree
        params = [{'name': 'Simulation Parameters', 'type': 'group', 'children': [
            {'name': 'birthrate','type':'float','value':self.world.birthrate},
            {'name': 'deathrate','type':'float','value':self.world.deathrate},
            {'name': 'regenerationrate','type':'float','value':self.world.regenerationrate},
            {'name': 'burdenrate','type':'float','value':self.world.burdenrate},
            {'name': 'economyaim','type':'float','value':self.world.economyaim},
            {'name': 'growthrate','type':'float','value':self.world.growthrate},
        ]}]
        self.params = Parameter.create(name='params', type='group', children=params)
        self.params.sigTreeStateChanged.connect(self.change)
        self.tree = ParameterTree()
        self.tree.setParameters(self.params, showTop=False)
        # show the app and the parameter tree
        self.tree.show()
        # add a slider
        #slider = win.addTickSliderItem()

        # start the main loop
        QtGui.QApplication.instance().exec_()

    def change(self,param, changes):
        print('##############')
        for parama, change, data in changes:
            path = self.params.childPath(parama)
            if data <= 0:
                self.logger.warn('{0} smaller 0 is not allowed'.format(path[-1]))
                setattr(self,path[-1],1)
            else:
                setattr(self,path[-1],data)
        self._r = r = self.world.solve()
        self.population_plot.setData(r['time'],r['population'])
        self.burden_plot.setData(r['time'],r['burden'])
        self.economy_plot.setData(r['time'],r['economy'])


if __name__ == '__main__':
    #main_old()
    app = pg.mkQApp()
    m = WorldSimpleGui
    app.exec_()
