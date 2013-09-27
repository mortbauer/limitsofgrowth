#coding=utf8
# virtualenv test_py2
import pandas
import numpy as np
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import logging
from scipy.integrate import ode, odeint
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
#from matplotlib import pyplot as plt


# state space representation

# y: population, burden, economy
# x: birth, death, ecocide, regeneration,  economicgrowth
# u: birthrate, deathrate, regenerationrate, burdenrate, economyaim, growthrate

class WorldSimpel(object):
    def __init__(self,birthrate=0.03,deathrate=0.01,regenerationrate=0.1,
                 burdenrate=0.02,economyaim=1,growthrate=0.05):
        self.birthrate = birthrate
        self.deathrate = deathrate
        self.regenerationrate = regenerationrate
        self.burdenrate = burdenrate
        self.economyaim = economyaim
        self.growthrate = growthrate
        self.time = np.linspace(0,100,1e5)
        # set up the logger
        self.logger = logging.getLogger(self.__class__.__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(levelname)s] %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        # set up the app and widgets
        win = pg.GraphicsWindow(title="very simple World Model")
        layout = QtGui.QGridLayout()
        win.setLayout(layout)
        text = QtGui.QLineEdit('enter text')
        layout.addWidget(text, 0, 1)
        plot = pg.PlotWidget()
        layout.addWidget(plot, 0, 0)
        #plot = win.addPlot(row=0,col=0)
        #plot.addLegend()
        # initial solution and plotting
        self._r = r = self.solve()
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
            {'name': 'birthrate','type':'float','value':self.birthrate},
            {'name': 'deathrate','type':'float','value':self.deathrate},
            {'name': 'regenerationrate','type':'float','value':self.regenerationrate},
            {'name': 'burdenrate','type':'float','value':self.burdenrate},
            {'name': 'economyaim','type':'float','value':self.economyaim},
            {'name': 'growthrate','type':'float','value':self.growthrate},
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

    def dx(self,x,time):
        population,burden,economy = x
        quality = burden**(-1)
        birth = self.birthrate * population * quality * economy
        death = population * self.deathrate * burden
        ecocide = self.burdenrate * economy * population
        regeneration = self.regenerationrate * burden if quality > 1 else self.regenerationrate
        economicgrowth = self.growthrate * economy * burden * (1-(economy*burden)/self.economyaim)
        return np.array([birth-death,ecocide-regeneration,economicgrowth])

    def solve(self,x0=[1,1,1]):
        res = pandas.DataFrame(odeint(self.dx,x0,self.time),
            columns=['population','burden','economy'])
        res['time'] = self.time
        return res

    def change(self,param, changes):
        for parama, change, data in changes:
            path = self.params.childPath(parama)
            if data <= 0:
                self.logger.warn('{0} smaller 0 is not allowed'.format(path[-1]))
                setattr(self,path[-1],1)
            else:
                setattr(self,path[-1],data)
        self._r = r = self.solve()
        self.population_plot.setData(r['time'],r['population'])
        self.burden_plot.setData(r['time'],r['burden'])
        self.economy_plot.setData(r['time'],r['economy'])

def f(x,t,birthrate=0.03,deathrate=0.01,
      regenerationrate=0.1,burdenrate=0.02,economyaim=1,growthrate=0.05):
    population,burden,economy = x
    quality = burden**(-1)
    birth = birthrate * population * quality * economy
    death = population * deathrate * burden
    ecocide = burdenrate * economy * population
    regeneration = regenerationrate * burden if quality > 1 else regenerationrate
    economicgrowth = growthrate * economy * burden * (1-(economy*burden)/economyaim)

    return np.array([birth-death,ecocide-regeneration,economicgrowth])

def solve():
    sol = []
    time = []
    solver = ode(lambda t,x:f(x,t))
    solver.set_initial_value([1,1,1],0)
    while solver.t < 500:
        solver.integrate(250,step=True)
        sol.append(solver.y)
        time.append(solver.t)
    return pandas.DataFrame(
        np.vstack(sol),index=pandas.Index(time,name='time'),
        columns=['population','burden','economy'])

def solve2(*args):
    t = np.linspace(0,100,1e5)
    res = pandas.DataFrame(
        odeint(f,[1,1,1],t,args=args),
        columns=['population','burden','economy'])
    res['time'] = t
    return res

def main_old():
    r = solve2()
    # matplotlib plot
    #sol.plot()
    # pyqtgraph plot
    app = pg.mkQApp()
    win = pg.GraphicsWindow(title="very simple World Model")
    layout = QtGui.QGridLayout()
    win.setLayout(layout)
    plot1 = pg.PlotWidget()#win.addPlot(title="Multiple curves")
    plot1.addLegend()
    plot1.plot(r['time'],r['population'],pen=(255,0,0),name='population')
    plot1.plot(r['time'],r['burden'],pen=(0,255,0),name='burden')
    plot1.plot(r['time'],r['economy'],pen=(0,0,255),name='economy')
    # text
    layout.addWidget(plot1, 0, 1, 3, 1)
    params = [{'name': 'Simulation Parameters', 'type': 'group', 'children': [
        {'name': 'economy','type':'float','value':10}]}]
    p = Parameter.create(name='params', type='group', children=params)
    def change(param, changes):
        print("tree changes:")
        for parama, change, data in changes:
            path = p.childPath(parama)
            if path is not None:
                childName = '.'.join(path)
            else:
                childName = parama.name()
            print('  parameter: %s'% childName)
            print('  change:    %s'% change)
            print('  data:      %s'% str(data))
            print('  ----------')

    p.sigTreeStateChanged.connect(change)
    t = ParameterTree()
    t.setParameters(p, showTop=False)
    t.show()
    app.exec_()


if __name__ == '__main__':
    #main_old()
    #app = pg.mkQApp()
    m = WorldSimpel()
    #app.exec_()
