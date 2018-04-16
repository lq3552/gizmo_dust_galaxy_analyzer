#!/usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import dust_gal_analyzer.DustyGalaxy as dg

def schechter(M,phis,Ms,alpha):
	x=M/Ms
	phi = phis*x**alpha*np.exp(-x)
	return np.log10(M),np.log10(phi)

def plot_obs(fname,fmt = '--o',color = 'k',label=''):
	data = np.log10(np.loadtxt(fname))
	x = data[:,0]
	y = data[:,1]
	yerr_l = y - data[:,2]
	yerr_h = data[:,3] - y
	plt.errorbar(x,y,yerr=[yerr_l,yerr_h],fmt=fmt,color = color,label=label)


num = ['003','006','010','013']#,'010']

plt.figure(figsize=(2*10,2*(8.0)))
for i in range(len(num)):
	gal_m25n256l = dg(glist=np.load('m25n256l/gal_snapshot_'+num[i]+'.npz'))
	gal_m25n256h = dg(glist=np.load('m25n256h/gal_snapshot_'+num[i]+'.npz'))
	if i < 3:
		gal_m25n256lld = dg(glist=np.load('m25n256lld/gal_snapshot_'+num[i]+'.npz'))
		gal_m25n256lhd = dg(glist=np.load('m25n256lhd/gal_snapshot_'+num[i]+'.npz'))
		gal_m25n256dh = dg(glist=np.load('m25n256dh/gal_snapshot_'+num[i]+'.npz'))
	gal_m25n256hld = dg(glist=np.load('m25n256hld/gal_snapshot_'+num[i]+'.npz'))
	gal_m25n256hhd = dg(glist=np.load('m25n256hhd/gal_snapshot_'+num[i]+'.npz'))

	ax=plt.subplot(2,2,i+1)
	plt.text(8.7,-2.5,'z ='+''+str(round(gal_m25n256l.get('redshift'),2)))
	
	l_m25n256l = gal_m25n256l.plot_dmf(sty = 'bo',label = '',xlabel='',ylabel='')
	l_m25n256h = gal_m25n256h.plot_dmf(sty = 'ro',label = '',xlabel='',ylabel='')
	if i < 3:
		l_m25n256lld = gal_m25n256lld.plot_dmf(sty = 'c^',label = '',xlabel='',ylabel='')
		l_m25n256lhd = gal_m25n256lhd.plot_dmf(sty = 'c*',label = '',xlabel='',ylabel='')
		l_m25n256dh = gal_m25n256dh.plot_dmf(sty = 'go',label = '',xlabel='',ylabel='')
	l_m25n256hld = gal_m25n256hld.plot_dmf(sty = 'm^',label = '',xlabel='',ylabel='')
	l_m25n256hhd = gal_m25n256hhd.plot_dmf(sty = 'm*',label = '',xlabel='',ylabel='')

	# observation
	if i == 1:
		plot_obs('dmf_z0.txt',label = 'Dunne+ 2003')
	if i == 2:
		plot_obs('dmf_z1.txt',label = 'Eales+ 2009')
	if i == 3:
		plot_obs('dmf_z1.txt',label = 'Dunne+ 2011')

	plt.axis([4.1,9.6,-5,-2])


	if i in [2,3]:
		plt.xlabel(r'$\log(M_d/M_\odot)$')
	else:
		ax.set_xticklabels([])
	if i in [0,2]:
		plt.ylabel(r'$\log\phi[{\rm Mpc^{-3}dex^{-1}}]$')
	else:
		ax.set_yticklabels([])

	if i == 0:
		legend2 = plt.legend([l_m25n256h,l_m25n256hld,l_m25n256hhd],['median shock','low shock','high shock'],loc=3)
		ax = plt.gca().add_artist(legend2)
		plt.legend([l_m25n256l,l_m25n256h,l_m25n256dh],['m25n256 low growth','m25n256 high growth','m25n256d high growth'],loc=4)
	if i == 1 :
		plt.legend(loc=4)

plt.subplots_adjust(left = 0.085, right = 0.99, bottom = 0.070, top = 0.98,hspace=0.0,wspace=0.0)
plt.savefig('dmf.png',dpi=300)



'''
# observation
Mdr= 10**np.linspace(6.1,8.7,20)
#m,phi = schechter(Md,2.97e-3,16.0e7,-1.01)
#plt.plot(m,phi,'o',label='z=0.2')
m,phi = schechter(Mdr,8.9e-4,47.0e7,-1.08) #doubtful
plt.plot(m,phi,'r--',label='z=2.0')
'''
