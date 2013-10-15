Simple World Model
##################

A simplified World Model for the Course `Kontinuierliche Simulation`, free
after J. W. Forrester's Model from 1971.

The package includes the `model.py` module which serves for the mathematical
stuff and some additional modules. `qt_gui.py` includes the most advanced gui
heavily utilizing `pyqtgraph`_ and `pyqt`_ and a second gui implemetation in
`kivy`_ which hasn't so many feautures as the qt one, but is also nice and
seems somehow more fun to program with. (This is anyways my first gui
programming so it was kind of a test.)

Installation
************
I suggest to setup a virtualenvironment, then one can do::

    pip install -e git+https://github.com/mortbauer/limitsofgrowth

One also needs quite many requirements, so also run::

    pip install -r requirements.txt

Depending on which gui version one wanna use you will also need::

    pip install -r requirements_kivy.txt

or::

    pip install -r requirements_qt.txt


GUI Usage
*********
To use the gui just run (in the activated virtualenv)::

    limitsofgrowth-qt

for the qt version, or::

    limitsofgrowth-kivy

for the kivy version.

Math Usage
**********
One can also just use the mathematic model without the gui. You basically do::

    from limitsofgrowth.model import WorldSimple

    w = WorldSimple()

To calculate the states now you can run::

    x = w.x_odeint(w.initial)

which will use the initial parameters for simulation.
