import dust_gal_analyzer.DustyGalaxy as dg
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.figsize'] = (10.0,8.0)
mpl.rcParams['lines.linewidth'] = 3.0
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['xtick.major.width'] = 2.0
mpl.rcParams['xtick.minor.width'] = 2.0
mpl.rcParams['ytick.major.width'] = 2.0
mpl.rcParams['ytick.minor.width'] = 2.0
mpl.rcParams['legend.fontsize'] = 20
mpl.rc('font',size=20)
import matplotlib.pyplot as plt

def SFRD_MD(z):
	return 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)

pa = 'm25n256l/'
pb = 'm25n256h/'

na = ['000','001','002','003','004','005','006','007','008','009','010','011','012','013']
nb = ['000','001','002','003','004','005','006','007','008','009','010','011','012','013']

za = []
zb = []
z = np.linspace(0,6,70)

SFRDa = []
SFRDb = []
rho_da = []
rho_db = []
rho_ga = []
rho_gb = []

for i in na:
	glist = np.load(pa+'gal_snapshot_'+i+'.npz')
	gal = dg(glist)
	SFRDa.append(gal.get('SFRD'))
	rho_da.append(gal._Od)
	rho_ga.append(gal._Og)
	za.append(gal.get('redshift'))

for i in nb:
	glist = np.load(pb+'gal_snapshot_'+i+'.npz')
	gal = dg(glist)
	SFRDb.append(gal.get('SFRD'))
	rho_db.append(gal._Od)
	rho_gb.append(gal._Og)
	zb.append(gal.get('redshift'))

za = np.array(za)
zb = np.array(zb)
SFRDa = np.array(SFRDa)
SFRDb = np.array(SFRDb)
rho_da = np.array(rho_da)
rho_db = np.array(rho_db)
rho_ga = np.array(rho_ga)
rho_gb = np.array(rho_gb)

plt.figure()
plt.plot(np.log10(1+za),np.log10(SFRDa),'ro',label='m25n256l')
plt.plot(np.log10(1+zb),np.log10(SFRDb),'bo',label='m25n256h')
plt.plot(np.log10(1+z),np.log10(SFRD_MD(z)),'k--',label='Maudau & Dickinson 2014')
plt.xlabel(r'$\log\ (1+z)$')
plt.ylabel(r'$\log$ SFRD $(M_\odot\ {\rm yr^{-1}\ Mpc^{-3}})$')
plt.legend()
plt.savefig('SFRDz.png',dpi=300)


plt.figure()
l1, = plt.plot(np.log10(1+za),np.log10(rho_da),'ro',label='m25n256l')
plt.plot(np.log10(1+zb),np.log10(rho_db),'bo',label='m25n256h')

l2, = plt.plot(np.log10(1+za),np.log10(rho_ga),'r*')
plt.plot(np.log10(1+zb),np.log10(rho_gb),'b*')

legend2 = plt.legend([l2,l1],['gas','dust'],loc=3)
ax = plt.gca().add_artist(legend2)
plt.xlabel(r'$\log\ (1+z)$')
plt.ylabel(r'$\log$ $\rho_{\rm met}$ $(M_\odot\ Mpc^{-3}})$')
plt.legend(loc = 4)
plt.savefig('rhoz.png',dpi=300)
