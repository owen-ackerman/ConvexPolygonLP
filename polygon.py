from matplotlib.patches import Polygon
import matplotlib.pyplot as plt

from matplotlib.widgets import Button
import numpy as np
from math import sqrt
from scipy.optimize import linprog
from datetime import datetime

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        ax.scatter(self.xs,self.ys)
        self.line.figure.canvas.draw()
    
    def finished(self, event):
        global start
        start=datetime.now()
        self.xs.append(0.01)
        self.ys.append(0.01)
        self.line.set_data(self.xs, self.ys)
        ax.scatter(self.xs,self.ys)
        LineEquation(self.xs,self.ys)
        # this data here will be used to construct the Matrix of line equations
        self.line.figure.canvas.draw()

#c = np.array
def LineEquation(xcoords, ycoords):
    global lines
    lines = np.zeros((len(xcoords)-1, 2))
    for i in range(len(xcoords)-1):
        a = ( (ycoords[i+1]-ycoords[i]) / (xcoords[i+1]-xcoords[i]) ) 
        # slope of line with coordinates (x1,y1), (x2,y2)
        b = ycoords[i] - xcoords[i]*a 
        # y intercept of line with coordinates (x1,y1), (x2,y2)
        lines[i, :] = [a,b]
        #populate array of a,b values of line equations in the form y=ax+b
    print("PRINTING LINES:")
    print(lines)
    #LinePlotter()
    EquationSeparator()
    OrientedLinePlotter()
        
    
def LinePlotter():
    x = np.linspace(0,6,100) 
    for i in range(len(lines)):
        ax.plot(x, lines[i][0]*x + lines[i][1])


def EquationSeparator():
    global upperBounds
    global lowerBounds
    global separator
    for i in range(len(lines)-1): #O(#lines)
        if lines[i][0] < 0:
            if lines[i+1][0] > 0:
                separator = i+1
                break

    upperBounds = np.zeros((separator, 2))
    upperBounds = lines[:separator]
    lowerBounds = np.zeros(((len(lines)-separator), 2))
    lowerBounds = lines[separator:]

    LinearSolver()

    

def OrientedLinePlotter():
    x = np.linspace(0,6,100) 
    for i in range(len(upperBounds)):
        ax.plot(x, upperBounds[i][0]*x + upperBounds[i][1], color='red')
    
    for i in range(len(lowerBounds)):
        ax.plot(x, lowerBounds[i][0]*x + lowerBounds[i][1], color='blue')

def FormattedUpperInequality():
    A = np.ones((len(upperBounds), 3)) #inequality variable coeficients
    B = np.ones((len(upperBounds), 1)) #inequality value
    for i in range(len(upperBounds)):
        a = sqrt(upperBounds[i][0]*upperBounds[i][0]+1)
        A[i][1] = -(upperBounds[i][0])/a
        A[i][2] = 1/a
        B[i][0] = upperBounds[i][1]/a

    return A, B

def FormattedLowerInequality():
    A = np.ones((len(lowerBounds), 3)) #inequality variable coeficients
    B = np.ones((len(lowerBounds), 1)) #inequality value
    for i in range(len(lowerBounds)):
        a = sqrt(lowerBounds[i][0]*lowerBounds[i][0]+1)
        A[i][1] = (lowerBounds[i][0])/a
        A[i][2] = -1/a
        B[i][0] = -lowerBounds[i][1]/a
    
    return A, B

def LinearSolver():
    c = np.array([-1,0,0])
    print(type(c))
    uB = FormattedUpperInequality()
    lB = FormattedLowerInequality()
    A_ub = np.concatenate((uB[0],lB[0]), axis = 0)
    b_ub = np.concatenate((uB[1],lB[1]), axis = 0)
    print("Printing A_ub")
    print(A_ub)
    print("Printing b_ub")
    print(b_ub)
    r_bounds = (0, None)
    s1_bounds = (0, 6.0)
    s2_bounds = (0, 6.0)
    bounds = [r_bounds, s1_bounds, s2_bounds]
    start=datetime.now()
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='interior-point')
    print(result)
    cc = plt.Circle((result.x[1],result.x[2]),result.x[0], alpha=0.3,color = 'b')
    #cc = plt.Circle((4,4),1)
    ax.add_patch(cc)
    print("ran")
    print(datetime.now()-start)




fig = plt.figure()
ax = fig.add_subplot(111)

plt.ylim(0,6)
plt.xlim(0,6)
ax.set_title('Please build a Convex Polygon')
line, = ax.plot([0.01], [0.01])  # empty line

axnext = plt.axes([.8, 0.05, 0.1, 0.075])
linebuilder = LineBuilder(line)
finished = Button(axnext, 'Finish')
finished.on_clicked(linebuilder.finished)

plt.show()

