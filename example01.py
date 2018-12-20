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
scat.set_cmap('jet')
nframes = 8
frames = []
for i in range(0,nframes):
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
    steel = pf.Material(30E6,36E3,0.284)
    element1 = pf.Element(1,node1,node2, steel, vertical_tube)
    element2 = pf.Element(2,node2,node3, steel, horizontal_tube)
    element3 = pf.Element(3,node3,node4, steel, vertical_tube)
    
    frame = pf.Structure([element1,element2,element3], [node1,node2,node3,node4])
    frame.solve()
    frame.calc_VMstress()
    frames.append(frame)
    
def animate(i):
    if i<nframes:
        scat.set_offsets(np.c_[frames[i].x_values,frames[i].y_values])
        scat.set_array(np.asarray(frames[i].von_misses))
    else:
        scat.set_offsets(np.c_[frames[2*(nframes-1)-i].x_values,frames[2*(nframes-1)-i].y_values])
        scat.set_array(np.asarray(frames[2*(nframes-1)-i].von_misses))
    

ani = FuncAnimation(fig, animate, frames=2*(nframes-1), interval=200, repeat=True) 
plt.show()