#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:30:30 2018

@author: ptdrow
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import planeFrame as pf

#animation changing an element length

fig, ax = plt.subplots(figsize=(5, 5))
ax.set(xlim=(-25, 150), ylim=(0, 130))

scat = ax.scatter([], [],s=1 )

def animate(i):
    node1 = pf.Node(1,0,0)
    node1.fix()
    
    node2 = pf.Node(2,i*0.625*12-1.25*12,10*12)
    node2.setForces(fx=10000)
    
    node3 = pf.Node(3,-i*0.625*12+11.25*12,10*12)
    node3.setForces(m=5000)
    
    node4 = pf.Node(4,10*12,0)
    node4.fix()
    
    vertical_tube = pf.Circular_Tube(6.449, 6.197)
    horizontal_tube = pf.Circular_Tube(4.647, 4.291)
    
    element1 = pf.Element(1,node1,node2, 30E6, vertical_tube)
    element2 = pf.Element(2,node2,node3, 30E6, horizontal_tube)
    element3 = pf.Element(3,node3,node4, 30E6, vertical_tube)
    
    frame = pf.Structure([element1,element2,element3], [node1,node2,node3,node4])
    frame.solve()
    frame.calc_VMstress()
    scat.set_offsets(np.c_[frame.x_values,frame.y_values])
    scat.set_array(np.asarray(frame.von_misses))
    scat.set_cmap('jet')

ani = FuncAnimation(fig, animate, frames=8, interval=100, repeat=True) 
plt.show()