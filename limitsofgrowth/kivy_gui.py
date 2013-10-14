#coding=utf8
# virtualenv test_py2
import sys
import pandas
import numpy as np
import logging
from scipy.integrate import ode, odeint

from colors import ColorWheel
from kivy.app import App
from kivy.graphics import Color,Mesh,Line
from kivy.uix.slider import Slider
from kivy.uix.label import Label
from kivy.uix.widget import Widget
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.relativelayout import RelativeLayout
from kivy.uix.gridlayout import GridLayout
from kivy.uix.tabbedpanel import TabbedPanelItem
from kivy.properties import NumericProperty,StringProperty, DictProperty, ObjectProperty, BooleanProperty, ReferenceListProperty, AliasProperty,ListProperty
from kivy.garden.graph import Graph, MeshLinePlot
from kivy.factory import Factory
from kivy.config import Config
from kivy.animation import Animation
from kivy.logger import Logger
from limitsofgrowth.model import WorldSimple
from limitsofgrowth.colorcube import print_color,rgb
from functools import partial

Config.set('graphics', 'width', '800')
Config.set('graphics', 'height', '500')

class TabbedPlot(TabbedPanelItem):
    pass

class Legend(GridLayout):
    labels = DictProperty({})

    def add_label(self,label):
        self.add_widget(label)

class LegendLabel(Widget):
    color = ListProperty()
    text = StringProperty()

class Plot(Graph):
    plotdict = DictProperty({})
    legend = ObjectProperty()

    def __init__(self,**kwargs):
        super(Plot,self).__init__(**kwargs)
        self.colorwheel = ColorWheel(start=0.6)

    def plot_xy(self,x,y,label):
        try:
            if not label in self.plotdict:
                c =[c/255 for c in self.colorwheel.next().rgb]
                self.plotdict[label] = d = {
                    'plot':MeshLinePlot(),
                    'label':LegendLabel(text=label,color=c),
                }
                d['plot']._color.rgb = c
                self.add_plot(d['plot'])
                self.legend.add_label(d['label'])
            self.plotdict[label]['plot'].points = np.array([x,y],copy=False).T.tolist()
            self.y_ticks_major = round(float(self.ymax - self.ymin)/10,1)
        except Exception as e:
            print('######################',e)

    def plot_dataframe(self,dataframe):
        ymaxold = 0
        xmaxold = 0
        yminold = 10e6
        xminold = 10e6
        for c in dataframe:
            x = dataframe.index
            y = dataframe[c]
            ymax = int(y.max())+1
            if ymax > ymaxold:
                ymaxold = ymax
            xmax = int(x.max())
            if xmax > xmaxold:
                xmaxold = xmax
            ymin = int(y.min())-1
            if ymin < yminold:
                yminold = ymin
            xmin = int(x.min())
            if xmin < xminold:
                xminold = xmin
            self.plot_xy(dataframe.index,dataframe[c],c)
        self.ymax = ymaxold
        self.xmax = xmaxold
        self.ymin = yminold
        self.xmin = xminold


class Controller(Widget):
    name = StringProperty()
    value = NumericProperty()
    max = NumericProperty()
    min = NumericProperty()

    def _update_value(self, value):
        try:
            self.value = round(float(value),3)
        except Exception as e:
            self.ids._text.text = str(self.value)

class LimitsOfGrowthApp(App):
    show_controllers = BooleanProperty(True)
    controllers_y = NumericProperty()
    world = ObjectProperty()

    def __init__(self,**kwargs):
        super(LimitsOfGrowthApp,self).__init__(**kwargs)
        self.register_event_type('on_parameters')
        self._parameters = {}
        self._controllers = {}

    def toggle_controllers(self):
        if self.show_controllers:
            self.controllers_y = self.root.ids.controllers.y
            y = - self.root.ids.controllers.height
        else:
            y = self.controllers_y

        self.show_controllers = not self.show_controllers

        Animation(y=y, d=.3, t='out_quart').start(
                self.root.ids.controllers)

    def build(self):
        # add a controller for each parameter
        i = 1
        for name,v in WorldSimple.params.items():
            c = Controller(
                name=name,value=v['initial'],max=v['max'],min=v['min'],
                y=(self.root.ids.controllers.height-10)*(1-float(i)/6),
                size_hint=(1,0.1)
            )
            c._update_value(v['initial'])
            #setattr(self,name,c)
            c.bind(value=self.set_parameters)
            self.root.ids.controllers.add_widget(c)
            self.root.ids[name] = c
            self._controllers[name] = c
            self._parameters[name] = v['initial']
            i+=1
        # create the world model
        self.world = WorldSimple()
        self.dispatch('on_parameters')

    def set_parameters(self,inst,value):
        if self._parameters[inst.name] != value:
            self._parameters[inst.name] = value
            self.dispatch('on_parameters')

    def new_plot(self):
        self.root.ids.tabbedpanel.add_widget(TabbedPlot(text='new plot'))

    @property
    def parameters(self):
        return self._parameters

    def on_parameters(self):
        self.root.ids.plot.plot_dataframe(self.world.x_odeint(self.parameters))


def main():
    LimitsOfGrowthApp().run()
if __name__ == '__main__':
    main()
