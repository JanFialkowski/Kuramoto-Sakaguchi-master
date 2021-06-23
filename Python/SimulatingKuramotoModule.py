# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:07:51 2020

@author: KingofPi
"""

import numpy as np
from scipy.integrate import odeint

def derivs(state,t, alpha = 0., beta = np.pi, eps = 0.01, sigma = 1., omega = 0.):
	N = int(-0.5+np.sqrt(0.25+len(state)))
	Phis = state[:N]
	Kappas = state[N:].reshape(N,N)
	DelPhis = (Phis-Phis[:,np.newaxis]).T
	DKappa = -eps*(np.sin(DelPhis+beta)+Kappas)
	DelPhis += alpha
	DelPhis = Kappas*np.sin(DelPhis)
	DPhis = DelPhis.sum(axis=1)*sigma/N
	DPhis += omega
	derivatives = np.concatenate((DPhis,DKappa.flatten()))
	return derivatives

def IntegratorwithKappas(state0, t, runs, alpha = 0., beta = np.pi, eps = 0.01, sigma = 1., omega = 0.):
	states = odeint(derivs,state0,t, args = (alpha,beta,eps,sigma,omega))
	Phis = [[] for i in xrange(int(-0.5+np.sqrt(0.25+len(state0))))]
	Kappa = [[],[]]
	for i in xrange(1,runs-1):
		for k in xrange(len(Phis)):
			for j in xrange(len(states)):
				Phis[k].append(states[j][k])
		for j in xrange(len(states)):
			Kappa[0].append(states[j][2*len(Phis)-1])
			Kappa[1].append(states[j][-len(Phis)])
		states = odeint(derivs, states[-1], t, args = (alpha,beta,eps,sigma,omega))
	print states
	for i in xrange(len(Phis)):
		for j in xrange(len(states)):
			Phis[i].append(states[j][i])
	for j in xrange(len(states)):
		Kappa[0].append(states[j][2*len(Phis)-1])
		Kappa[1].append(states[j][-len(Phis)+1])
	return Phis,Kappa

def Integrator(state0, t, runs, alpha = 0., beta = np.pi, eps = 0.01, sigma = 1., omega = 0.):
	Phis = [[] for i in xrange(int(-0.5+np.sqrt(0.25+len(state0))))]
	for i in xrange(runs):
		states = odeint(derivs,state0,t, args = (alpha,beta,eps,sigma,omega))
		for k in xrange(len(Phis)):
			for j in xrange(len(states)):
				Phis[k].append(states[j][k])
		state0 = states[-1]
	return np.array(Phis)