import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.animation as manimation
import Kuramoto as K
from mpl_toolkits.mplot3d import Axes3D
#from scipy.integrate import solve_ivp

alpha = 0.4*np.pi
beta = -0.9*np.pi
eps = 0.01

def Gamma(n):
	delta = 2*(alpha+beta)
	s = np.sin(delta)
	c = np.cos(delta)
	G = n*s
	G2 =-((1+c)*n-1)
	Out=np.arctan2(G,G2)
	return Out

def derivs(state,t):
	N = int(-0.5+np.sqrt(0.25+len(state)))
	Phis = state[:N]
	Kappas = state[N:].reshape(N,N)
	DelPhis = (Phis-Phis[:,np.newaxis]).T
	DKappa = -eps*(np.sin(DelPhis+beta)+Kappas)
	DelPhis += alpha
	DelPhis = Kappas*np.sin(DelPhis)
	DPhis = DelPhis.sum(axis=1)/Phis.shape[0]*(-1)
	derivatives = np.concatenate((DPhis,DKappa.flatten()))
	if t%100==0:
		print t/100.
	return derivatives

def Calccoords(state):
	while state[0]>2*np.pi:
		state[0]-=2*np.pi
	while state[1]>2*np.pi:
		state[1]-=2*np.pi
	while state[2]>2*np.pi:
		state[2]-=2*np.pi
	state.sort()
	s0 = 0
	s1 = state[1]-state[0]
	s2 = state[2]-state[0]
	l3 = (2*np.pi-s2)/(2*np.pi)
	l2 = np.abs(s2-s1)/(2*np.pi)
	l1 = np.abs(s1-s0)/(2*np.pi)
	Sum = l1+l2+l3
	return [l1,l2,l3]

def Coordstopointinspace(l):
	x = (l[1]+l[0]*0.5)/(l[0]+l[1]+l[2])
	y = l[0]*np.sqrt(3)/2/(l[0]+l[1]+l[2])
	return [x,y]

def CalcR2(Phis):
	return np.sum((np.exp(2j*Phis)),axis=0)/len(Phis)

def Integrator(state0, t, runs):
	states = odeint(derivs,state0,t)
	Phis = [[] for i in xrange(int(-0.5+np.sqrt(0.25+len(state0))))]
	Kappa = [[],[]]
	for i in xrange(1,runs-1):
		for k in xrange(len(Phis)):
			for j in xrange(len(states)):
				Phis[k].append(states[j][k])
		for j in xrange(len(states)):
			Kappa[0].append(states[j][2*len(Phis)-1])
			Kappa[1].append(states[j][-len(Phis)])
		states = odeint(derivs,states[-1],t)
	print states
	for i in xrange(len(Phis)):
		for j in xrange(len(states)):
			Phis[i].append(states[j][i])
	for j in xrange(len(states)):
		Kappa[0].append(states[j][2*len(Phis)-1])
		Kappa[1].append(states[j][-len(Phis)])
	return Phis,Kappa

if __name__=="__main__":
	N = 3
	N2 = 1.
	Phi0 = np.zeros(N)#np.random.rand(50)*2*np.pi*0.01
	G = Gamma(N2/N)+alpha+beta
	Phi0[0]= -(G)
#	Phi0[1] = np.pi
	Kappa0 = -np.sin((Phi0-Phi0[:,np.newaxis]).T+beta).flatten()
	Phi0[0] += 0.1
	state0 = np.concatenate((Phi0,Kappa0))
	t = np.arange(0,30000,1)
	print (-G)/np.pi
	Phis, Kappa = Integrator(state0, t, 1)
	state0[0]-=0.2
	Phis2, Kappa2 = Integrator(state0, t, 1)
	state0[0]+=0.1
	#"""
	state0[2*len(Phi0)-2:2*len(Phi0)]+=0.01
	Phis3, Kappa3 = Integrator(state0, t, 1)
	state0[2*len(Phi0)-2:2*len(Phi0)]-=0.02
	Phis4, Kappa4 = Integrator(state0, t, 1)
	state0[2*len(Phi0)-2:2*len(Phi0)]-=0.01
	state0[-len(Phi0)]+=0.01
	state0[-2*len(Phi0)] +=0.01
	Phis5, Kappa5 = Integrator(state0, t, 1)
	state0[-len(Phi0)]-=0.02
	state0[-2*len(Phi0)] -=0.02
	Phis6, Kappa6 = Integrator(state0, t, 1)
	#"""
#begin the plotting
	plt.plot(Kappa[0])
	plt.plot(Kappa[1])
	plt.show()
	fig = plt.figure()
	ax = fig.add_subplot(111, projection = "3d")
	#ax2 = fig.add_subplot(122,projection = "3d")
	ax.plot(Kappa3[0], Kappa3[1], (np.array(Phis3[0])-np.array(Phis3[1]))/np.pi, label = "K12+")
	ax.plot(Kappa6[0], Kappa6[1], (np.array(Phis6[0])-np.array(Phis6[1]))/np.pi, label = "k21-")
	ax.plot(Kappa5[0], Kappa5[1], (np.array(Phis5[0])-np.array(Phis5[1]))/np.pi, label = "K21+")
	ax.plot(Kappa4[0], Kappa4[1], (np.array(Phis4[0])-np.array(Phis4[1]))/np.pi, label = "K12-")
	ax.plot(Kappa[0], Kappa[1], (np.array(Phis[0])-np.array(Phis[1]))/np.pi, label = "T+")
	ax.scatter(-np.sin(-G-Phi0[1]+np.pi+beta), -np.sin(Phi0[1]+G-np.pi+beta), (-G-Phi0[1])/np.pi-1, c="black")
	ax.scatter(-np.sin(-G-Phi0[1]+np.pi+beta), -np.sin(Phi0[1]+G-np.pi+beta), (-G-Phi0[1])/np.pi+1, c="black")
	ax.scatter(-np.sin(-G-Phi0[1]+beta), -np.sin(Phi0[1]+G+beta), (-G-Phi0[1])/np.pi, c="black")
	ax.scatter(-np.sin(beta),-np.sin(beta),0,c="orange")
	ax.scatter(-np.sin(beta+np.pi),-np.sin(beta+np.pi),1,c="orange")
	ax.scatter(-np.sin(beta),-np.sin(beta),2,c="orange")
	ax.plot(Kappa2[0], Kappa2[1], (np.array(Phis2[0])-np.array(Phis2[1]))/np.pi, label = "T-")
	ax.set_xlabel("Kappa12")
	ax.set_ylabel("Kappa21")
	ax.set_zlabel("Theta")
	ax.legend()
	#ax2.set_xlabel("Kappa12")
	#ax2.set_ylabel("Kappa21")
	#ax2.set_zlabel("Theta")
	plt.show()
	Phisets = np.array([Phis,Phis2,Phis3,Phis4,Phis5,Phis6])
	for thing in Phisets:
		plt.plot(np.mod(thing[:,-1],2*np.pi), label = np.where(Phisets == thing))
	plt.legend()
	plt.show()
	fig, axs = plt.subplots(ncols=2, sharey = True)
	ax = axs[0]
	ax2 = axs[1]
	ax.set_xlabel("Time step")
	ax2.set_xlabel("Time step")
	ax.set_title("Positive perturbation")
	ax2.set_title("Negative perturbation")
	ax.set_ylabel("Phase")
#	ax2.set_ylabel("phase")
	for i in xrange(len(Phis)):
		ax.plot(Phis[i])#Phasen-time-plot
		ax2.plot(Phis2[i])
	plt.show()
	APhis = np.mod(np.array(Phis),2*np.pi)
	APhis2 = np.mod(np.array(Phis2),2*np.pi)
	R = CalcR2(np.array(Phis))
	R2 = CalcR2(np.array(Phis2))
	print np.abs(R)
	print np.abs(R2)
	plt.subplot(121, polar=True, title = "Positive Perturbation").plot(np.angle(R),np.abs(R))#Orderparameter time plot
	plt.subplot(122, polar=True, title = "Negative Perturbation").plot(np.angle(R2),np.abs(R2))#Orderparameter time plot
	plt.show()
	#APhis.sort(axis=0)
	#APhis2.sort(axis=0)
	plt.plot(APhis.T[-1]/np.pi)#last Phases over index
	plt.plot(APhis2.T[-1]/np.pi)#last Phases over index
	plt.show()
	ax = plt.subplot(121, polar=True)
	ax2 = plt.subplot(122, polar=True)
	for i in xrange(len(APhis)):
		ax.scatter(np.angle(np.exp(1j*APhis[i][-1])),1)#Oscillators on circle
		ax2.scatter(np.angle(np.exp(1j*APhis2[i][-1])),1)#Oscillators on circle
	plt.show()
	fig = plt.figure()
	ax = fig.add_subplot(121)
	ax2 = fig.add_subplot(122)
	ax.plot((np.array(Phis[0])-np.array(Phis[1]))/np.pi, label = "Theta/pi")
	ax.plot(np.abs(R), label = "R")
	ax2.plot((np.array(Phis2[0])-np.array(Phis2[1]))/np.pi, label = "Theta/pi")
	ax2.plot(np.abs(R2), label = "R")
	plt.legend()
	plt.show()
	alpha = np.linspace(0,np.pi/2, num=1000)[:,np.newaxis]
	beta = (np.linspace(0,2*np.pi,num=4000)-np.pi)
	"""
	a,b = K.L1(1/50.,alpha,beta,eps)
	c,d = K.L2(1/50.,alpha,beta,eps)
	e,f = K.L3(1/50.,alpha,beta,eps)
	CountMatrix = 48*(a.real>=0)+48*(b.real>=0)+(e.real>=0)+(f.real>=0)
	plt.imshow(CountMatrix, origin = "lower", extent = (-1,1,0,0.5))
	plt.colorbar()
	plt.show()
	
	a,b = K.L1(1/3.,alpha,beta,eps)
	c,d = K.L2(1/3.,alpha,beta,eps)
	e,f = K.L3(1/3.,alpha,beta,eps)
	CountMatrix= 1*(a.real>=0)+1*(b.real>=0)+0*(c.real>=0)+0*(d.real>=0)+(e.real>=0)+(f.real>=0)
	plt.imshow(CountMatrix, origin = "lower", extent = (-1,1,0,0.5))
	plt.colorbar()
	plt.show()
	"""
"""
	FFMpegWriter = manimation.writers['ffmpeg']
	metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
	writer = FFMpegWriter(fps=25, metadata=metadata)

	fig = plt.figure()
	ax = fig.add_subplot(121)
	ax2 = fig.add_subplot(122, polar = True)
	for Set in Phis:
		ax.plot(Set)
	ax2.plot(np.angle(R),np.abs(R))
	with writer.saving(fig, "Testing.mp4", 100):
		for t in range(len(R)):
			ax.plot([t,t],[np.min(Phis),np.max(Phis)], "black")
			ax2.scatter(np.angle(R[t]),np.abs(R[t]),c="orange")
			writer.grab_frame()
			del ax.lines[-1]
			del ax2.collections[-1]
			print t
"""
