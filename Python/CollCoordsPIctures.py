# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:02:10 2020

@author: KingofPi
"""

import numpy as np
import matplotlib.pyplot as plt
import SimulatingKuramotoModule as S

N = 100
alpha= 0
beta = np.pi/2
eps = 0.01
sigma = 2.55


def Order(Phis):
	Dummies = 1.j*Phis
	return np.abs(np.exp(Dummies).sum(axis=0)/Phis.shape[0])

Phis = np.zeros(N)
Kappas = -np.sin(Phis[:,np.newaxis]-Phis+beta)
t = np.linspace(0,10,num = 10)
omega = np.linspace(-1,1,num=N)
state0 = np.concatenate((Phis,Kappas.flatten()))
TPhis = S.Integrator(state0,t,100,alpha = alpha, beta = beta, omega=omega, sigma = sigma, eps = eps)
print(TPhis)
r = Order(TPhis)
print(r)
plt.plot(r)
plt.show()
plt.plot(np.remainder(TPhis[:,r.argmax()], 2*np.pi))
plt.plot(np.remainder(TPhis[:,r.argmin()], 2*np.pi))
plt.show()