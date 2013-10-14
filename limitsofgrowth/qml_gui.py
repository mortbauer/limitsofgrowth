import os
import sys
from pkg_resources import resource_filename
from PyQt5 import QtCore, QtGui, QtQuick, QtWidgets

from limitsofgrowth.model import WorldSimple

class QWorldModel(QtCore.QObject):
    def __init__(self):
        super(QWorldModel,self).__init__()
        self._world = WorldSimple()

    @QtCore.pyqtProperty(list)
    def population(self):
        return QtQuickself._world.x_odeint(self._world.initial)['population'].values


# instantiate example
world = QWorldModel()

# run app
app = QtWidgets.QApplication(sys.argv)
view = QtQuick.QQuickView()
view.rootContext().setContextProperty('world', world)
view.setSource(QtCore.QUrl(resource_filename(__name__,'qml/limitsofgrowth.qml')))
view.show()
app.exec_()
