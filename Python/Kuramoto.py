import numpy as np
import matplotlib.pyplot as plt

def Gamma(n, alpha, beta):
	delta = 2*(alpha+beta)
	s = np.sin(delta)
	c = np.cos(delta)
	G = n*s
	G2 =-((1+c)*n-1)
	Out=np.arctan2(G,G2)
	return Out

def derivs(state,t, alpha, beta, eps):
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

def L1(n, alpha, beta, eps):
	Psi = -Gamma(n, alpha, beta)-alpha-beta
	Theta = alpha + beta
	Thate = alpha-beta
	p = 0.5*(np.sin(Thate)-(1-n)*np.sin(Theta)-n*np.sin(-2*Psi+Theta)+2*eps)
	q = -eps*((1-n)*np.sin(Theta)+n*np.sin(-2*Psi+Theta))
	p2 = p**2
	return 0.5*p+np.sqrt(0.25*p2-q+0j), 0.5*p-np.sqrt(0.25*p2-q+0j)

def L2(n, alpha, beta, eps):
	Psi = -Gamma(n, alpha, beta)-alpha-beta
	Theta = alpha + beta
	Thate = alpha-beta
	p = 0.5*(np.sin(Thate)-(1-n)*np.sin(2*Psi+Theta)-n*np.sin(Theta)+2*eps)
	q = -eps*((1-n)*np.sin(2*Psi+Theta)+n*np.sin(Theta))
	p2 = p**2
	return 0.5*p+np.sqrt(0.25*p2-q+0j), 0.5*p-np.sqrt(0.25*p2-q+0j)

def L3(n, alpha, beta, eps):
	Psi = -Gamma(n, alpha, beta)-alpha-beta
	Theta = alpha + beta
	Thate = alpha-beta
	p = 0.5*(np.sin(Thate)-(1-n)*np.sin(2*Psi+Theta)-n*np.sin(-2*Psi+Theta)+2*eps)
	q = -eps*((1-n)*np.sin(2*Psi+Theta)+n*np.sin(-2*Psi+Theta))
	p2 = p**2
	return 0.5*p+np.sqrt(0.25*p2-q+0j), 0.5*p-np.sqrt(0.25*p2-q+0j)
"""
def Jacobian(state, alpha, beta, eps):
	N = -0.5+np.sqrt(0.25-len(state))
	Phis = state[:N]
	Kappas = state[N:].reshape(N,N)
	DelPhis = Phis[:,np.newaxis]-Phis
	Jacob[:N,:N] = 1/N*Kappas*np.cos(DelPhis+alpha)
	for i in xrange(N):
		J[i,i] -= J.sum(axis=1)
	
	Jacob = np.zeros((len(state),len(state)))
	Jacob[N:,N:]=np.eye(N**2)*(-eps)
"""
if __name__ == "__main__":
	alpha = np.linspace(0,np.pi/2, num=1000)[:,np.newaxis]
	beta = (np.linspace(0,2*np.pi,num=4000)-np.pi)
	eps = 0.01
	n= 1./50.
	L1=L1(n,alpha, beta, eps)
	L2=L2(n,alpha, beta, eps)
	L3=L3(n,alpha, beta, eps)
	fig = plt.figure()
	ax = fig.add_subplot(321)
	im = ax.imshow(L1[0].real==0.,origin = "lower")
	ax2 = fig.add_subplot(322)
	im2 = ax2.imshow(L1[1].real==0.,origin = "lower")
	ax3 = fig.add_subplot(323)
	im3 = ax3.imshow(L2[0].real==0.,origin = "lower")
	ax4 = fig.add_subplot(324)
	im4 = ax4.imshow(L2[1].real==0.,origin = "lower")
	ax5 = fig.add_subplot(325)
	im5 = ax5.imshow(L3[0].real==0.,origin = "lower")
	ax6 = fig.add_subplot(326)
	im6 = ax6.imshow(L3[1].real==0.,origin = "lower")
	fig.colorbar(im, ax=ax)
	fig.colorbar(im2, ax=ax2)
	fig.colorbar(im3, ax=ax3)
	fig.colorbar(im4, ax=ax4)
	fig.colorbar(im5, ax = ax5)
	fig.colorbar(im6, ax = ax6)
	plt.show()
