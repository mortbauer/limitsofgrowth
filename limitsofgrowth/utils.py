import pandas
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

class ScaledDataFrame(pandas.DataFrame):
    def __init__(self,*args ,**kwargs):
        super(ScaledDataFrame,self).__init__(*args,**kwargs)
        self._scale()

    @property
    def scales(self):
        return self._scales

    def _scale(self):
        if not hasattr(self,'_scales'):
            self._scales = {c:10**round(log(self[c].max(),10)) for c in self.columns}
            for c in self.columns:
                self[c] = self[c]/self._scales[c]

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

