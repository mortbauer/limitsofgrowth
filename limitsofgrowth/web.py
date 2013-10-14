from limitsofgrowth import WorldSimple
from pygal import Line
from flask import Flask
import numpy as np

app = Flask(__name__)
app.debug = True
w = WorldSimple()


@app.route("/simulation.svg")
def simulation():
    x = w.x_odeint(w.initial)
    line = Line(show_dots=False,show_minor_x_labels=False)
    for c in x:
        line.add(c,x[c])
    line.x_labels = map(str,x.index.values)
    line.x_labels_major = map(str,range(x.index[0],x.index[-1]+50,50))
    return line.render_response()

if __name__ == "__main__":
    app.run()
