#coding=utf8
# virtualenv test_py2
import sys
import time

from colors import ColorWheel
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

import logging
import gettext
from pkg_resources import resource_filename

from limitsofgrowth.model import WorldSimple

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

    def plot_all(self,dataframe,over=None,postfix='',prefix=''):
        if over:
            index = dataframe[over]
        else:
            index = dataframe.index
        for column in dataframe:
            if column != over:
                self.plot_single(
                    index,dataframe[column],'%s%s%s'%(prefix,column,postfix))

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

class WorldSimpleGui(QtGui.QMainWindow):
    """ segementation fault bug, see:
        https://groups.google.com/d/msg/pyqtgraph/juoqAiXABYY/BYpvCQeWSgMJ
    """
    def __init__(self):
        super(WorldSimpleGui,self).__init__()
        # absolutly main
        self.dataframeplots = {}
        self.data_tree = {}
        self.world = WorldSimple(resolution=1000)
        self.parameters = {param:value['initial']
                            for param,value in self.world.params.items()}
        self.dataframeplots['Realtime Simulation'] = DataFramePlot()
        self.dataframeplots['Realtime Simulation'].plot_all(
            self.world.x_odeint(self.parameters))
        # set main window gui
        self.setCentralWidget(self.dataframeplots['Realtime Simulation'].plotwidget)
        self.setWindowTitle("Simple Limits of Growth World Model")
        self.mainplot = self.dataframeplots['Realtime Simulation']
        self.resize(1000,600)
        self.tools = QtGui.QTabWidget()
        dock = QtGui.QDockWidget('Controllers')
        dock.setWidget(self.tools)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea,dock)
        # additional
        self.create_parameter_widget()
        self.create_main_controls_widget()
        self.create_data_widget()
        self._add_plot_to_data_widget('Realtime Simulation')
        self.opt_basinhopping = [  1.00227541e-02,   9.36331825e+00,   1.00000000e+01,
         2.23513971e-02,   9.93493979e+01,   5.95912808e-01]

    def _remove_plot_to_data_widget(self,name):
        self.data_widget.removeTopLevelItem(self.data_tree.pop(name)['item'])

    def _add_plot_to_data_widget(self,name):
        def make_callback(plot,curve):
            return lambda q_point: self.data_tree_context_menu(q_point,plot,curve)
        plot = self.dataframeplots[name]
        if not name in self.data_tree:
            tree_item = pg.TreeWidgetItem()
            tree_label = QtGui.QLabel(name)
            tree_label.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
            tree_label.customContextMenuRequested.connect(make_callback(name,None))
            tree_item.setWidget(0,tree_label)
            self.data_tree[name] = {'item':tree_item,'children':{}}
            self.data_widget.addTopLevelItem(tree_item)
        else:
            tree_item = self.data_tree[name]['item']
        for curve_name,curve in plot.plots.items():
            if not curve_name in self.data_tree[name]['children']:
                child = pg.TreeWidgetItem()
                label = QtGui.QLabel(curve_name)
                label.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
                label.customContextMenuRequested.connect(make_callback(name,curve_name))
                tree_item.addChild(child)
                self.data_widget.setItemWidget(child,0,label)
                self.data_tree[name]['children'][curve_name] = child

    def data_tree_context_menu(self,q_point,plotname,curve):
        menu = QtGui.QMenu(self)
        delete_action = menu.addAction('clear')
        plot = self.dataframeplots[plotname]
        legend = plot.plotwidget.plotItem.legend
        def clear():
            if not curve:
                plot.plotwidget.clear()
                # https://groups.google.com/d/msg/pyqtgraph/UTwLwC5mQnQ/HxGgXc0xzbMJ
                legend.items = []
                legend.updateSize()
                plot.plots = {}
                self._remove_plot_to_data_widget(plotname)
            else:
                plot.plotwidget.removeItem(plot.plots.pop(curve))
                legend.items = [x for x in legend.items if not x[1].text == curve]
                legend.updateSize()
                self.data_tree[plotname]['item'].removeChild(self.data_tree[plotname]['children'].pop(curve))

        delete_action.triggered.connect(clear)
        move_action = QtGui.QWidgetAction(menu)
        move = QtGui.QWidget()
        move_button = QtGui.QPushButton('move to')
        layout = QtGui.QHBoxLayout()
        plot_selector = QtGui.QComboBox()
        plot_selector.setEditable(True)
        for item in self.dataframeplots:
            plot_selector.addItem(item)
        layout.addWidget(move_button)
        layout.addWidget(plot_selector)
        move.setLayout(layout)
        move_action.setDefaultWidget(move)
        menu.addAction(move_action)
        menu.exec_(self.data_widget.mapToGlobal(q_point))


    def create_data_widget(self):
        self.data_widget = pg.TreeWidget()
        self.data_widget.resize(400,self.data_widget.height())
        self.data_widget.setColumnCount(1)
        self.data_widget.setHeaderLabels(['plots'])
        self.tools.addTab(self.data_widget,'Data Explorer')

    def create_real_data_tree_widget(self):
        self.real_dw = pg.DataTreeWidget(data={'optimization':self.optimized})

    def create_main_controls_widget(self):
        widget = QtGui.QWidget()
        widget.resize(400,widget.height())
        layout = QtGui.QVBoxLayout()
        sens_group = QtGui.QGroupBox('Sensitivity Analysis')
        sens_group_layout = QtGui.QVBoxLayout()
        sens_button = QtGui.QPushButton('scipy.odeint simultaneous')
        sens_button.clicked.connect(self.plot_sensitivity)
        sens_button_cv = QtGui.QPushButton('assimulo.cvode simultaneous')
        sens_button_cv.clicked.connect(self.plot_sensitivity_cvode)
        sens_button_cvn = QtGui.QPushButton('assimulo.cvode natural')
        sens_button_cvn.clicked.connect(self.plot_sensitivity_cvode_natural)
        ref_group = QtGui.QGroupBox('Reference Data')
        ref_layout = QtGui.QVBoxLayout()
        ref_button = QtGui.QPushButton('world indicators worldbank')
        ref_button.clicked.connect(lambda : self.plot_dataframe(
            self.world.reference_data/self.world.reference_data.loc[1960],
            postfix=' ref'))
        ref_layout.addWidget(ref_button)
        opt_button = QtGui.QPushButton('plot optimized')
        opt_button.clicked.connect(self.plot_optimized)
        ref_layout.addWidget(opt_button)
        diff_button = QtGui.QPushButton('plot difference')
        diff_button.clicked.connect(self.plot_diff)
        ref_layout.addWidget(diff_button)
        real_group = QtGui.QGroupBox('Simulation')
        real_group_layout = QtGui.QVBoxLayout()
        self.birth_mod = QtGui.QCheckBox('use modified birth')
        real_group_layout.addWidget(self.birth_mod)
        real_group.setLayout(real_group_layout)
        sens_group_layout.addWidget(sens_button)
        sens_group_layout.addWidget(sens_button_cv)
        sens_group_layout.addWidget(sens_button_cvn)
        sens_group.setLayout(sens_group_layout)
        ref_group.setLayout(ref_layout)
        # select the plot
        self.plot_selector = QtGui.QComboBox(self)
        self.plot_selector.setEditable(True)
        for item in self.dataframeplots:
            self.plot_selector.addItem(item)
        layout.addWidget(self.plot_selector)
        layout.addWidget(sens_group)
        layout.addWidget(ref_group)
        layout.addWidget(real_group)
        widget.setLayout(layout)
        self.tools.addTab(widget,'Plot Creation')

    def create_parameter_widget(self):
        self.sliders = {}
        widget = QtGui.QWidget()
        widget.resize(400,widget.height())
        layout = QtGui.QVBoxLayout()
        for param in self.world.params:
            value = self.world.params[param]
            slider = Parameter(param,value=value['initial'],
                               xmin=value['min'],xmax=value['max'],step=0.01)
            slider.valueChanged.connect(self.change)
            layout.addWidget(slider)
            self.sliders[param] = slider
        default_button = QtGui.QPushButton('Reset Parameters')
        default_button.clicked.connect(self.reset_parameters)
        layout.addWidget(default_button)
        optimized_button = QtGui.QPushButton('Set to Optimized')
        optimized_button.clicked.connect(self.set_optimized)
        layout.addWidget(optimized_button)
        widget.setLayout(layout)
        self.tools.addTab(widget,'Parameters')

    def reset_parameters(self):
        for param,value in self.world.params.items():
            self.sliders[param].value = value['initial']

    @property
    def optimized(self):
        if not hasattr(self,'_optimized'):
            QtGui.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
            self._optimized = (self.opt_basinhopping,)#self.world.fit_with_data_bfgs()
            QtGui.QApplication.restoreOverrideCursor()
        return self._optimized

    def plot_diff(self):
        self.plot_dataframe(self.world.diff(self.world.x_odeint(self.optimized[0])),postfix=' diff')

    def plot_optimized(self):
        x = self.world.x_odeint(self.optimized[0])
        x = x/x.loc[1960]
        self.plot_dataframe(x.loc[self.world.reference_data.index],postfix=' optimized')

    def set_optimized(self):
        for i,param in enumerate(self.world.params):
            self.sliders[param].value = self.optimized[0][i]

    def plot_dataframe(self,dataframe,postfix='',prefix=''):
        selected_plot_name = str(self.plot_selector.currentText())
        if not selected_plot_name in self.dataframeplots:
            self.dataframeplots[selected_plot_name] = DataFramePlot(
                window_title=selected_plot_name)
            self.dataframeplots[selected_plot_name].plotwidget.closeEvent = lambda event:self._remove_plot_to_data_widget(selected_plot_name)
        self.dataframeplots[selected_plot_name].plot_all(dataframe,postfix=postfix,prefix=prefix)
        self.dataframeplots[selected_plot_name].plotwidget.show()
        self._add_plot_to_data_widget(selected_plot_name)

    def plot_sensitivity(self):
        try:
            s = self.world.s_odeint(self.parameters)
            data = s[['{0},economyaim'.format(x) for x in self.world.cols]]
            self.plot_dataframe(data,postfix=' odeint')
        except:
            logger.warn('no solution for sensitivities found')

    def plot_sensitivity_cvode_natural(self):
        try:
            s = self.world.s_cvode_natural(self.parameters)
            data = s[['{0},economyaim'.format(x) for x in self.world.cols]]
            self.plot_dataframe(data,postfix=' cvode natural')
        except:
            logger.warn('no solution for sensitivities found')

    def plot_sensitivity_cvode(self):
        try:
            s = self.world.s_cvode(self.parameters)
            data = s[['{0},economyaim'.format(x) for x in self.world.cols]]
            self.plot_dataframe(data,postfix=' cvode')
        except:
            logger.warn('no solution for sensitivities found')

    def matplotlib_plot(self):
        # not in use currently
        self.p = MatplotlibPlot()
        self.p.dataframe = data
        self.p.on_draw()
        self.p.show()

    def change(self,param, value):
        oldvalue = self.parameters[str(param)]
        self.parameters[str(param)] = value
        t1 = time.time()
        inpv = self.world.create_input_vector(self.parameters)
        if self.birth_mod.isChecked():
            func=self.world.create_dx(inpv,birth_callback=self.world.birth_mod)
        else:
            func=self.world.create_dx(inpv)
        res = self.world.x_odeint(func=func)
        if res:
            self.plot_dataframe(res)
            logger.info('recalculated in {0}s'.format(time.time()-t1))
        else:
            logger.warn('failed with parameter "{0}" equal {1}'.format(param,value))


def main():
    m = WorldSimpleGui()
    m.show()
    qt_app.exec_()

def get_reference_data():
    # no data available
    water_pollution = ['EE.BOD.%s.ZS'% sector for sector in
                       ('CGLS','FOOD','MTAL','OTHR','PAPR','TXTL','WOOD')]

    indicators = {'population':'SP.POP.TOTL',
                  'economy':'NY.GDP.MKTP.CD',
                  'burden':'EN.ATM.CO2E.KT'}

    countries = {'world':'WLD',
                 'high income':'HIC',
                 'low income':'LIC',
                 'low & middle income':'LMY',
                 'high middel income':'HIC',
                 'low middle income':'LMC'}

    wc = WorldBankClient()

    for country in countries.values():
        for indicator in indicators:
            if type(indicators[indicator]) == list:
                for ind in indicators[indicator]:
                    wc.get_indicator_by_country(ind,country)
            else:
                wc.get_indicator_by_country(indicators[indicator],country)
    return wc


if __name__ == '__main__':
    main()
