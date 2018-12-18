#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 11:32:34 2018

@author: Pedro Villarroel @ptdrow
"""
import math
import numpy as np
import matplotlib.pyplot as plt

class Node:
    count = 0
    nodes = []
    def __init__(self, number, x_coordinate, y_coordinate):
        self.number = number
        self.x = x_coordinate
        self.y = y_coordinate
        Node.nodes.append(self)
        Node.count += 1
        self.setForces()
        self.setDofs()
    
    def setForces(self, fx=0, fy=0, m=0):
        self.forces = np.array([fx, fy, m])
    
    def setDofs(self, x = True, y = True, angle = True):
        self.dofs = [x, y, angle]
    
    def fix(self):
        self.setDofs(False,False,False)


class Cross_Section:
    def __init__(self, area, moment_of_inertia):
        self.A = area
        self.I = moment_of_inertia
    
class Circular_Tube(Cross_Section):
    def __init__(self, radius_external, radius_internal):
        self.re = radius_external
        self.ri = radius_internal
        self.A = math.pi * (pow(self.re,2) - pow(self.ri,2))
        self.I = math.pi * (pow(self.re,4) - pow(self.ri,4)) / 4
    
    def calc_Qdivb(self,y):
        y1 = abs(y)
        if y1 <= self.ri:
            Q = 2/3*(math.sqrt(pow(self.re**2 - y1**2 , 3)) - 
                     math.sqrt(pow(self.ri**2 - y1**2 , 3)))
            b = 2*(math.sqrt(self.re**2 - y1**2 ) - 
                     math.sqrt(self.ri**2 - y1**2 ))
            return Q/b
        else:
            Qdivb = 1/3 * (self.re**2 - y1**2 ) #for avoiding errors at y1 = re because b=0
            return Qdivb

    
class Element:
    count = 0
    elements = []
    
    def __init__(self, number, first_node, second_node,
                 Young_module, cross_section):
        self.number = number
        self.E = Young_module
        self.cSection = cross_section
        self.setNodes(first_node, second_node)
        self.calcVariables()
        self.createKLocal()
        Element.elements.append(self)
        Element.count += 1
        self.displacements = np.zeros(6)
        self.forces = np.zeros(6)
    
    def calcVariables(self):
        self.delta_x = self.node2.x - self.node1.x
        self.delta_y = self.node2.y - self.node1.y
        self.L = math.sqrt(pow(self.delta_x , 2) + pow(self.delta_y , 2))
        self.alfa = math.atan2(self.delta_y , self.delta_x)
        self.c1 = self.cSection.A * self.E / self.L
        self.c2 = self.E * self.cSection.I / pow(self.L,3)
    
    def setNodes(self, new_node1, new_node2):
        self.node1 = new_node1
        self.node2 = new_node2
    
    def createKLocal(self):
        k1 = self.c1
        k2 = self.c2 * 12
        k3 = self.c2 * self.L * 6
        k4 = self.c2 * pow(self.L , 2) * 4
        k5 = self.c2 * pow(self.L , 2) * 2
        
        # Local stiffness matrix in local coordinates 
        self.KLocal_e = np.asarray([[ k1, 0 , 0 ,-k1, 0 , 0 ], # 1st row
                                      [ 0 , k2, k3, 0 ,-k2, k3], # 2nd row
                                      [ 0 , k3, k4, 0 ,-k3, k5], # 3rd row
                                      [-k1, 0 , 0 , k1, 0 , 0 ], # 4th row
                                      [ 0 ,-k2,-k3, 0 , k2,-k3], # 5th row
                                      [ 0 , k3, k5, 0 ,-k3, k4], # 6th row
                                      ])
        C = math.cos(self.alfa)
        S = math.sin(self.alfa)
        self.T = np.asarray([[ C, S, 0, 0, 0, 0], # 1st row
                        [-S, C, 0, 0, 0, 0], # 2nd row
                        [ 0, 0, 1, 0, 0, 0], # 3rd row
                        [ 0, 0, 0, C, S, 0], # 4th row
                        [ 0, 0, 0,-S, C, 0], # 5th row
                        [ 0, 0, 0, 0, 0, 1], # 6th row
                        ])
        # Local stiffness matrix in global coordinates
        self.KLocal = np.linalg.multi_dot([np.transpose(self.T), self.KLocal_e, self.T])
    
    def moment(self, x):
        m1 = self.forces[2]
        m2 = self.forces[5]
        M = (m2 + m1)/self.L * x - m1
        return M
    
    def normal_stress(self, x, y):
        M = self.moment(x)
        fx2 = self.forces[3]
        A = self.cSection.A
        I = self.cSection.I
        sigmaX = fx2 / A - M * y /I # positive M gives positive stress (traction) with negative y
        return sigmaX
    
    def shear_stress(self, y):
        V = self.forces[4]
        return(V * self.cSection.calc_Qdivb(y) / self.cSection.I)
        
    def von_misses_stress(self, x, y):
        sigmaX = self.normal_stress(x, y)
        taoXY = self.shear_stress(y)
        return(math.sqrt(sigmaX**2 / 2 + 3 * taoXY**2))
        
        
class Structure:
    def __init__(self, elements, nodes):
        self.elements = elements
        self.nodes = nodes
        self.e_count = len(elements)
        self.n_count = len(nodes)
        
    def solve(self):
        self.assembly_KGlobal()
        self.apply_BC()
        self.solve_displacements()
        self.solve_local_forces()
        
    def assembly_KGlobal(self):
        self.KGlobal = np.zeros((self.n_count * 3, self.n_count * 3)) # 3 dofs for every node

        #Assembly KLocales in KGlobal
        for e in self.elements:
            for i in range(0, 6):
                if i < 3:
                    iGlobal = (e.node1.number - 1) * 3 + i
                else:
                    iGlobal = (e.node2.number - 2) * 3 + i
        
                for j in range(0, 6):
                    if j < 3:
                        jGlobal = (e.node1.number - 1) * 3 + j
                    else:
                        jGlobal = (e.node2.number - 2) * 3 + j
            
                    self.KGlobal[iGlobal, jGlobal] += e.KLocal[i,j]
        
    def apply_BC(self):
        self.dofs = self.nodes[0].dofs
        self.forces = self.nodes[0].forces
        for node in self.nodes[1:]:
            self.dofs = self.dofs + node.dofs
            self.forces = np.hstack((self.forces, node.forces))
            
    def solve_displacements(self):
        KReduced = self.KGlobal[self.dofs][...,self.dofs]
        FReduced = self.forces[self.dofs]
        self.displacements = np.zeros(len(self.dofs))
        self.displacements[self.dofs] = np.dot(np.linalg.inv(KReduced),FReduced)
        
    def solve_local_forces(self):
        for e in Element.elements:
            for i in range (0, 3):
                e.displacements[ i ]   = self.displacements[(e.node1.number-1)*3 + i] 
                e.displacements[3 + i] = self.displacements[(e.node2.number-1)*3 + i]
            e.forces = np.linalg.multi_dot([e.KLocal_e, e.T, e.displacements])
            
    def calc_VMstress(self, n=18,r=10):
        m = n*r
        von_misses = np.zeros((n+1)*(m+1)*self.e_count)
        x_values = np.zeros((n+1)*(m+1)*self.e_count)
        y_values = np.zeros((n+1)*(m+1)*self.e_count)
        k = -1
        for e in self.elements:
            u = - e.L/m
            for i in range(0, m+1):
                u += e.L/m
                v = (n+2)/n*e.cSection.re
                for j in range(0, n+1):
                    k +=1
                    v -= 2*e.cSection.re/n
                    x_values[k] = e.node1.x + u * e.delta_x/e.L - v * e.delta_y/e.L
                    y_values[k] = e.node1.y + u * e.delta_y/e.L + v * e.delta_x/e.L
                    von_misses[k] = e.von_misses_stress(u, v)        
        
        all_values = list(zip(x_values,y_values,von_misses))
        all_values = sorted(all_values, key=lambda values: values[2])
        self.x_values, self.y_values, self.von_misses = zip(*all_values)
        
    def plot_VMstress(self):
        plt.clf()
        plt.scatter(self.x_values, self.y_values, c=self.von_misses, 
                    cmap='jet', s=1 )
        plt.title('Von Misses stress for the structure')
        plt.ylabel('y')
        plt.xlabel('x')
        plt.colorbar()
        plt.show()

