#coding=utf8
# virtualenv test_py2
import os
import sys
import time
from math import log
import re
import glob
import json
import traversedata
import requests

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
from scipy import optimize
from numpy import linalg
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

    params_list = ['birthrate','deathrate','regenerationrate',
                            'burdenrate','economyaim','growthrate']
    cols = ['population','burden','economy']

    reference_data_names = {'population':'SP.POP.TOTL',
                            'economy':'NY.GDP.MKTP.CD',
                            'burden':'EN.ATM.CO2E.KT'}

    def __init__(self, resolution=250):
        self.time = np.arange(1900,1900+250+1,1)
        self.sensitivities = []
        self.initial = {param:v['initial'] for param,v in self.params.items()}
        for column in self.cols:
            for param in self.params:
                self.sensitivities.append('{0},{1}'.format(column,param))

    def create_dx(self,*args):
        params = args[0]
        if type(params)== dict:
            birthrate = params['birthrate']
            deathrate = params['deathrate']
            regenerationrate = params['regenerationrate']
            burdenrate = params['burdenrate']
            economyaim = params['economyaim']
            growthrate = params['growthrate']
        elif len(params) == 6:
            birthrate = params[0]
            deathrate = params[1]
            regenerationrate = params[2]
            burdenrate = params[3]
            economyaim = params[4]
            growthrate = params[5]
        else:
            raise Exception('wrong number of arguments')
        def func(time,x):
            """
            the change of the system variables population, burden and economy
            x: [population,burden,economy]
            """
            population,burden,economy = x
            quality = x[1]**(-1)
            birth = birthrate * population * quality * economy
            death = population * deathrate * burden
            ecocide = burdenrate * economy * population
            if quality > 1:
                regeneration = regenerationrate * burden
            else:
                regeneration = regenerationrate
            economicgrowth = growthrate * economy * burden \
                    * (1-(economy*burden)/economyaim)
            f = lambda : [birth-death,ecocide-regeneration,economicgrowth]
            fnachx = lambda : np.array([
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
            fnachp = lambda : np.array([
                [population/burden*economy,
                - population+burden,0,0,0,0],
                [0,0,-burden if burden<1 else -1,
                economy*population,0,0],
                [0,0,0,0,growthrate*economy**2*burden**2/economyaim**2,
                economyaim*burden*(1-(economy*burden/economyaim))]
            ])
            return {'f':f,'fx':fnachx,'fp':fnachp}
        return func

    def x_odeint(self,params):
        func = self.create_dx(params)
        res,info = odeint(lambda x,t:func(t,x)['f'](),[1.0,1.0,1.0], self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe = pandas.DataFrame(res,
                columns=['population','burden','economy'],index=self.time)
            return dataframe

    def x_cvode(self,params):
        func = self.create_dx(params)
        problem = Explicit_Problem(lambda t,x:func(t,x)['f'](), [1.0,1.0,1.0],0)
        sim = CVode(problem)
        t,x = sim.simulate(250,len(self.time)-1)
        dataframe = pandas.DataFrame(x,
                columns=['population','burden','economy'],index=self.time)
        return dataframe

    def create_ds(self,*args):
        f = self.create_dx(*args)
        dx = np.empty((21,))
        def func(t,x):
            _f = f(t,x[:3])
            dx[:3] = _f['f']()
            _s = x[3:].reshape((3,6))
            dx[3:] = (_f['fx']().dot(_s)+_f['fp']()).reshape((18,))
            return dx
        return func

    def s_odeint(self,params):
        func = self.create_ds(params)
        s0 = np.ones(21)
        s0[3:]=0
        res,info = odeint(lambda s,t:func(t,s),s0, self.time,
                          full_output=True,printmessg=False,mxhnil=0)
        if info['message'] == "Integration successful.":
            dataframe = pandas.DataFrame(res,columns=self.cols+self.sensitivities)
            return dataframe

    def s_cvode(self,params):
        func = self.create_ds(params)
        s0 = np.ones(21)
        s0[3:]=0
        problem = Explicit_Problem(lambda t,s:func(t,s),s0,1900)
        sim = CVode(problem)
        t,s = sim.simulate(1900+250,self.time.shape[0]-1)
        dataframe = pandas.DataFrame(
            s,columns=self.cols+self.sensitivities,index=self.time)
        return dataframe

    def s_cvode_natural(self,params):
        problem = Explicit_Problem(lambda t,x,p:self.create_dx(p)(t,x)['f'](),
                                   [1.0,1.0,1.0],0,
                                   [params[p] for p in self.params_list])
        sim = CVode(problem)
        sim.report_continuously = True
        t,x = sim.simulate(250,self.time.shape[0]-1)
        dataframe = pandas.DataFrame(x,
                columns=['population','burden','economy'])
        d = {}
        sens = np.array(sim.p_sol)
        for i,col in enumerate(self.cols):
            for j,param in enumerate(
                ('birthrate','deathrate','regenerationrate',
                 'burdenrate','economyaim','growthrate')):
                d['{0},{1}'.format(col,param)] = sens[j,:,i]
        dataframe_sens = pandas.DataFrame(d,index=self.time)
        return dataframe_sens

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
        sol = pandas.DataFrame(np.empty((t1/dt,22)),columns=self.cols+sensitivities)
        i = 0
        #return solver_x,solver_s
        while solver_x.successful() and solver_s.successful() and solver_x.t < t1:
            solver_x.integrate(solver_x.t+dt)
            sol.iloc[i][self.cols] = solver_x.y
            solver_s.set_f_params(params,solver_x.y)
            solver_s.integrate(solver_s.t+dt)
            sol.iloc[i][sensitivities] = solver_s.y
            i += 1
        return sol

    @property
    def reference_data(self):
        if not hasattr(self,'_reference_data'):
            wc = WorldBankClient()
            wc.load_from_hdf(resource_filename(__name__,'data/world_indicators.hdf'))
            d = wc.indicators_by_countries('WLD')
            self._reference_data = d.rename_axis(
                {v:key for key,v in self.reference_data_names.items()})
        return self._reference_data

    def diff(self,result):
        return result.loc[self.reference_data.index]/result.loc[1960]-\
                self.reference_data/self.reference_data.loc[1960]

    def create_residum_func(self):
        initial = [self.initial[x] for x in self.params_list]
        bnds = [(0,None) for i in range(len(self.params))]
        def residum(p):
            try:
                result = self.x_odeint(p)
                diff = self.diff(result)
                return linalg.norm(diff)
            except:
                return 10e6
        return residum,initial,bnds

    def fit_with_data(self):
        func,initial,bnds = self.create_residum_func()
        r = optimize.minimize(func,initial,method='SLSQP',
                              bounds=bnds,options={'maxiter':10e4,'ftol':10})
        return r


class WorldBankClient(object):
    _REGION_CODES = 'http://worldbank.270a.info/classification/region.html'
    BASE_URL = 'http://api.worldbank.org/'
    _ID = re.compile(r'world_bank_data_(?P<country>[A-Z]*)_(?P<indicator>[A-Z.0-9]*).json')

    def __init__(self):
        self.indicators = {}
        self._plots = {}

    def indicators_by_countries(self,country):
        indicators = self.indicators.keys()
        d = {ind:v[country] for ind,v in self.indicators.items() if country in v}
        return pandas.DataFrame(d)

    @property
    def all_countries(self):
        if not hasattr(self,'_all_countries'):
            self._all_countries = {x['name']:x['id'] for x in self.countries_basicinfo}
        return self._all_countries

    @property
    def countries_basicinfo(self):
        if not hasattr(self,'_countries_basicinfo'):
            self._countries_basicinfo = self._get(
                '{0}countries/all/'.format(self.BASE_URL))
        return self._countries_basicinfo

    @property
    def indicators_basicinfo(self):
        if not hasattr(self,'_indicators_basicinfo'):
            self._indicators_basicinfo = self._get(
                '{0}indicators/all/'.format(self.BASE_URL))
        return self._indicators_basicinfo

    @property
    def all_indicators(self):
        if not hasattr(self,'_all_indicators'):
            self._all_indicators = {x['name']:x['id'] for x in self.indicators_basicinfo}
        return self._all_indicators

    def get_indicator_by_country(self,indicator,country):
        if not indicator in self.indicators or not country in self.indicators[indicator]:
            d = traversedata.WorldBankData(
                self._get_indicator_by_country(indicator,country),
                country=country,indicator=indicator)
            if indicator in self.indicators:
                self.indicators[indicator][country] = d.data
            else:
                self.indicators[indicator] = pandas.DataFrame(d.data,columns=[country])
        return self.indicators[indicator][country]

    def load_from_json_and_dir(self,path):
        for filepath in glob.glob('{0}world_bank_data_*.json'.format(path)):
            with open(filepath,'r') as f:
                m = self._ID.match(f.name.split('/')[-1])
                indicator = m.group('indicator')
                country = m.group('country')
                if not indicator in self.indicators or not country in self.indicators[indicator]:
                    d = traversedata.WorldBankData(
                        json.loads(f.read()))
                    if indicator in self.indicators:
                        self.indicators[indicator][country] = d
                    else:
                        self.indicators[indicator] = pandas.DataFrame(d.data,columns=[country])

    def save_to_json(self,path):
        with open(path,'w') as f:
            f.write(json.dumps({ind:v.raw for ind,v in self.indicators.items()}))

    def save_to_hdf(self,path,mode='a'):
        store = pandas.HDFStore(path,mode=mode)
        try:
            for indicator,countries in self.indicators.items():
                for country,d in countries.items():
                    d.data.to_hdf(store,'{0}/{1}'.format(
                        indicator.replace('.','_'),country))
        finally:
            store.close()

    def load_from_hdf(self,path):
        store = pandas.HDFStore(path,mode='r')
        try:
            for key in store.keys():
                indicator,country= key.lstrip('/').split('/')
                indicator = indicator.replace('_','.')
                if indicator in self.indicators:
                    self.indicators[indicator][country] = store[key]
                else:
                    self.indicators[indicator] = pandas.DataFrame(store[key],columns=[country])
        finally:
            store.close()


    def _get_indicator_by_country(self,indicator,country):
        return self._get('{0}countries/{1}/indicators/{2}'.format(
            self.BASE_URL,country,indicator))

    def _get(self,url,page=1,per_page=50,until=None):
        result = []
        def get(page):
            r = requests.get(
                '{0}?format=json&per_page={1}&page={2}'.format(url,per_page,page))
            if r.ok:
                status,res = r.json()
                result.extend(res)
                if page < status['pages'] and until and page < until:
                    get(page+1)
        get(1)
        return result

    def plot(self):
        self._plots['main'] = DataFramePlot()
        self._plots['main']

class DataFramePlot(object):
    def __init__(self,window_title='',size=(1000,600)):
        self.plots = {}
        self.colorwheel = ColorWheel(start=0.6)
        self.plotwidget = pg.PlotWidget()
        self.plotwidget.addLegend()
        self.plotwidget.setWindowTitle(window_title)
        self.plotwidget.resize(*size)
        self.plotwidget.show()

    def plot_single(self,x,y,name,label=None):
        if not label:
            label = name
        if name in self.plots:
            self.plots[name].setData({'x':x,'y':y},name=_(label))
        else:
            self.plots[name] = self.plotwidget.plotItem.plot({'x':x,'y':y},
                pen=tuple(self.colorwheel.next().rgb),name=_(label)
            )

    def plot_all(self,dataframe,over=None):
        if over:
            index = dataframe[over]
        else:
            index = dataframe.index
        for column in dataframe:
            if column != over:
                self.plot_single(index,dataframe[column],column)

class ScaledDataFrame(pandas.DataFrame):
    def __init__(self,*args,**kwargs):
        super(ScaledDataFrame,self).__init__(*args,**kwargs)
        self._scale()

    @property
    def scales(self):
        return self._scales

    def _scale(self):
        if not hasattr(self,'_scales'):
            self._scales = {c:10**round(log(self[c].max(),10)) for c in self.columns}
            for c in self.columns:
                self[c]/=self._scales[c]

    def plot(self,widget=None,over=None):
        if not widget:
            widget = DataFramePlot()
        if over:
            index = self[over]
        else:
            index = self.index.values.astype(np.float)
        for column in self:
            if column != over:
                widget.plot_single(index,self[column],column,label='{0} * {1:.1e}'
                                   .format(column,1.0/self.scales[column]))
        return widget

class Parameter(QtGui.QGroupBox):
    def __init__(self,name,value=1,xmin=0,xmax=10,step=1.0):
        super(Parameter,self).__init__(_(name))
        self.param_name = name
        self._value = value
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

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self,value):
        if value != self._value:
            self._value = value
            self.textbox.setText(str(value))
            self.slider.setValue((value-self.min)/self.step)
            self.valueChanged.emit(self.param_name,value)

    def on_value_slider_changed(self,value):
        self.value = self.step*value+self.min

    def on_value_textbox_changed(self):
        value = str(self.textbox.text())
        if value:
            self.value = float(value)

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

class WorldSimpleGui(QtGui.QMainWindow):
    """ segementation fault bug, see:
        https://groups.google.com/d/msg/pyqtgraph/juoqAiXABYY/BYpvCQeWSgMJ
    """
    def __init__(self):
        super(WorldSimpleGui,self).__init__()
        # absolutly main
        self.dataframeplots = {}
        self.world = WorldSimple(resolution=1000)
        self.parameters = {param:value['initial']
                            for param,value in self.world.params.items()}
        self.dataframeplots['main'] = DataFramePlot()
        self.dataframeplots['main'].plot_all(
            self.world.x_odeint(self.parameters))
        # set main window gui
        self.setCentralWidget(self.dataframeplots['main'].plotwidget)
        self.setWindowTitle("Simple Limits of Growth World Model")
        self.resize(1000,600)
        self.tools = QtGui.QTabWidget()
        dock = QtGui.QDockWidget('Controllers')
        dock.setWidget(self.tools)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea,dock)
        # additional
        self.create_parameter_widget()
        self.create_main_controls_widget()

    def create_main_controls_widget(self):
        widget = QtGui.QWidget()
        widget.resize(300,widget.height())
        layout = QtGui.QVBoxLayout()
        sens_button = QtGui.QPushButton('Sensitivity Analysis')
        sens_button.clicked.connect(self.plot_sensitivity)
        sens_button_cvn = QtGui.QPushButton('Sensitivity Analysis Cvode Natural')
        sens_button_cvn.clicked.connect(self.plot_sensitivity_cvode_natural)
        default_button = QtGui.QPushButton('Reset Parameters')
        default_button.clicked.connect(self.reset_parameters)
        layout.addWidget(sens_button)
        layout.addWidget(sens_button_cvn)
        layout.addWidget(default_button)
        widget.setLayout(layout)
        self.tools.addTab(widget,'Tasks')

    def create_parameter_widget(self):
        self.sliders = {}
        widget = QtGui.QWidget()
        widget.resize(300,widget.height())
        layout = QtGui.QVBoxLayout()
        for param,value in self.world.params.items():
            slider = Parameter(param,value=value['initial'],
                               xmin=value['min'],xmax=value['max'],step=0.01)
            slider.valueChanged.connect(self.change)
            layout.addWidget(slider)
            self.sliders[param] = slider
        widget.setLayout(layout)
        self.tools.addTab(widget,'Parameters')

    def reset_parameters(self):
        for param,value in self.world.params.items():
            self.sliders[param].value = value['initial']

    def plot_sensitivity(self):
        try:
            s = self.world.s_odeint(self.parameters)
            data = s[['{0},economyaim'.format(x) for x in self.world.cols]]
            if not 'sensitivities' in self.dataframeplots:
                self.dataframeplots['sensitivities'] = DataFramePlot(window_title='Sensitivity Analysis')
            self.dataframeplots['sensitivities'].plot_all(data)
            self.dataframeplots['sensitivities'].plotwidget.show()
        except:
            logger.warn('no solution for sensitivities found')
    def plot_sensitivity_cvode_natural(self):
        try:
            s = self.world.s_cvode_natural(self.parameters)
            data = s[['{0},economyaim'.format(x) for x in self.world.cols]]
            if not 'sensitivities_cvode_natural' in self.dataframeplots:
                self.dataframeplots['sensitivities_cvode_natural'] = DataFramePlot(window_title='Sensitivity Analysis Cvode Natural')
            self.dataframeplots['sensitivities_cvode_natural'].plot_all(data)
            self.dataframeplots['sensitivities_cvode_natural'].plotwidget.show()
        except:
            logger.warn('no solution for sensitivities found')

        #self.p = MatplotlibPlot()
        #self.p.dataframe = data
        #self.p.on_draw()
        #self.p.show()

    def change(self,param, value):
        oldvalue = self.parameters[str(param)]
        self.parameters[str(param)] = value
        t1 = time.time()
        res = self.world.x_odeint(self.parameters)
        if res:
            self.dataframeplots['main'].plot_all(res)
            logger.info('recalculated in {0}s'.format(time.time()-t1))
        else:
            logger.warn('failed with parameter "{0}" equal {1}'.format(param,value))


def main():
    m = WorldSimpleGui()
    m.show()
    qt_app.exec_()

if __name__ == '__main__':
    main()
