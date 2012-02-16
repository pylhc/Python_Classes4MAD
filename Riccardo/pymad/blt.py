#!/bin/sh

""":"
exec python $0 ${1+"$@"}
"""

# ---------------------- HelloWavesExt.py ----------------------
# 
# This program demonstrates how to make an animation of the 
# wave equation. The user can specify the initial 
# conditions by dragging the graph.
# 
# This is a slightly extended version of the HelloWaves.py script,
# also containing a short description of the numerical model used
# to simulate the motion of a guitar string.

from Tkinter import *        # The Tk package
import Pmw                   # The Python MegaWidget package
import math                  # import the sin-function

master = Tk()                # build Tk-environment
npoints = 15                 # use 10 points on each curve

def mouseDrag(event):
    global idx

    u[idx] = g.yaxis_invtransform(event.y)
    g.element_configure("string", ydata = tuple(u))

def mouseDown(event):
    global idx

    g.element_bind("string", "<Motion>", mouseDrag)
    el = g.element_closest(event.x, event.y)
    idx  = el["index"]

def mouseUp(event):
    g.element_unbind("string", "<Motion>")
    init()

def step():
    global u, up, um, C

    # Description of the mathematical model used for simulating
    # waves on a guitar string:
    #
    # Waves on a string are described (as many other phenomena
    # like light and radio waves) by the wave equation
    #
    #            u_tt = u_xx
    #
    # i.e. the second derivative in time of a function u(x,t) equals
    # the second derivative in space (sometimes you see the wave velocity
    # in here as well, but we have now scaled it away).
    #
    # The physical interpretation of u in the current context is the 
    # displacement of the string from its equilibrium position.
    #
    # The equation u_tt = u_xx is to be solved for 0<x<1.5 in this example.
    # At time t=0 we need to specify the shape of the string (this can
    # be done interactively!) and the velocity of the string, which is
    # naturally set to zero.
    #
    # The mathematical model outlined so far is then discretized by
    # a numerical method, here the finite difference method. This means
    # that the second-order derivatives are replaced by finite differences
    # according to the formula:
    #
    #      f''(x[i]) = (f(x[i+1]) -2f(x[i]) +f(x[i-1]))/(dx*dx)
    #
    # where dx is the distance between x[i-1] and x[i], and x[i] and x[i+1].
    # If we imagine that all f values are stored in an array f[i], we can
    # write the approximation like this:
    #
    #      f''(x[i]) = (f[i+1] -2*f[i] + f[i-1])/(dx*dx)
    #
    # To apply such a formula to our present equation, we need to store
    # the computed values of u in an array, which for simplicity can be
    # written as u[i,j], where i is a counter for points in the x direction
    # and j is a counter for points in t direction.
    # Inserting the approximations to the second derivatives in the
    # equation u_tt = u_xx results in a recurrence relation
    #
    # (u[i,j+1] -2*u[i,j] + u[i,j-1])/(dt*dt) =
    #               (u[i+1,j] -2*u[i,j] + u[i-1,j])/(dx*dx)
    #
    # This formula can be solved with respect to u[i,j+1], which 
    # is the new value of u we compute (assuming that u at time levels
    # j and j-1 are already computed).
    #
    # u[i,j+1] = 2*u[i,j] - u[i,j-1] + C*C*(u[i+1,j] -2*u[i,j] + u[i-1,j])
    #
    # with C*C = dt*dt/(dx*dx).
    #
    # It is only necessary to store three time levels: j+1, j, and j-1,
    # i.e., all previous levels can be thrown away. To this end, we introduce
    # three arrays u[i], up[i], and um[i], for u[i,j], u[i,j+1], and
    # u[i,j-1]. Moreover, we need a special value of u[i,-1] to get
    # started with the scheme (this value incorporates the start condition
    # that the velocity of the string is zero):
    #
    #       um[i] = u[i] + 0.5*C*C*(u[i+1] - 2*u[i] + u[i-1])
    #
    # The parameter C is important: C>1 will introduce numerical instability
    # and the displacement of the string will explode (try it!).
    # C=1 is a very special case where the numerical model (which is
    # usually only an approximation to the mathematical model) is exact
    # for any choice of dx=dt. However, if you make a wavy initial condition,
    # the exact solution obtained by C=1 is less attractive from a
    # visual point of view, so we have set C=0.8 here to introduce some
    # numerical errors that make the graph look nicer ;-)
    #
    # As u=0 always at the end points, we only update up (the new values)
    # at the interior x points 1...(npoints-1)
    
    # here is the magic line that finds new displacements of the
    # points along the string:
    
    for i in range(1, len(u)-1):
        up[i] = 2*u[i] -um[i] + C*C*(u[i+1] - 2*u[i] + u[i-1])
        
    # shuffle arrays (to be ready for calling step() at the next time level):
    um = u[:]
    u = up[:]

    # update the graph:
    g.element_configure("string", ydata = tuple(u))

    
def run():
    g.element_configure("string", symbol="")
    
    ntimesteps = 200  # the length of the animation
    for t in range(ntimesteps):
       step()
       master.after(20)           # wait 0.02 second
       master.update_idletasks()  # update screen

       g.element_configure("string", symbol="circle")


def init():
    global u, up, um, C
    C = 0.8    
    # make special value for um (incorporating du/dt=0 at t=0):
    um = u[:]  # ensure that um is an array of correct length...
    for i in range(1, len(um)-1):
        um[i] = u[i] + 0.5*C*C*(u[i+1] - 2*u[i] + u[i-1])

    up = u[:]

def zero():
    for i in range(len(u)):
        u[i] = 0
        g.element_configure("string", ydata = tuple(u))
        init()

if not Pmw.Blt.haveblt(master):     # Is Blt installed?
    print("BLT is not installed!")

else:
    vector_x = []   # x coordinates of the grid points
    u        = []   # string displacement at the x values
    
    # make a default graph
    for x in range(npoints+1):
        vector_x.append(x*0.1)
        if(x < 2.0*npoints/3):
            u.append(x/(2.0*npoints/3))
        else:
            u.append(float(npoints-x)/(npoints/3))

    init()
    g = Pmw.Blt.Graph(master)                  
    g.pack(expand=1, fill='both')

    curvename = 'string'                       
    g.line_create(curvename,                
                  xdata=tuple(vector_x), 
                  ydata=tuple(u),        
                  pixels=7,
                  smooth='natural')      # smooth the curve by splines

    g.element_bind("string", "<ButtonPress>",   mouseDown)
    g.element_bind("string", "<ButtonRelease>", mouseUp  )

    g.yaxis_configure(min=-1, max=1)
    g.configure(title='Hello world of waves', width=250, height=200)

    buttons = Pmw.ButtonBox(master, labelpos='n', label_text='Options')
    buttons.pack(fill='x', expand=0, padx=10, pady=10)

    buttons.add('Simulate',   command=run)
    buttons.add('Step',       command=step)
    buttons.add('Zero',       command=zero)
    buttons.add('Quit',       command=master.quit)
    
    master.mainloop()
   

