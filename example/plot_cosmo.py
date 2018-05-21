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
import sys

def SFRD_MD(z):
	return 0.015*(1+z)**2.7/(1+((1+z)/2.9)**5.6)

p = sys.argv[1] + '/'
n = ['000','001','002','003','004','005','006','007','008','009','010','011','012']#,'013']
#n = ['000']
z = []
zt = np.linspace(0,6,70)

SFRD = []
rho_g = []
grho_g = []

for i in n:
	glist = np.load(p+'gal_snapshot_'+i+'.npz')
#	glist = np.load('gal_snap_m25n256_135.npz')
#	glist = np.load('m25n256lch/gal_snapshot_013.npz')
	gal = dg(glist)
	SFRD.append(gal.get('SFRD'))
	rho_g.append(gal.get('gas_metal_density'))
	grho_g.append(np.sum(gal.get('gas_mass')*gal.get('gas_Z'))/gal.get('volume'))
	z.append(gal.get('redshift'))

print grho_g
z = np.array(z)
SFRD = np.array(SFRD)
rho_g = np.array(rho_g)

plt.figure()
plt.plot(np.log10(1+z),np.log10(SFRD),'ro',label=p[:-1])
plt.plot(np.log10(1+zt),np.log10(SFRD_MD(zt)),'k--',label='Maudau & Dickinson 2014')
plt.xlabel(r'$\log\ (1+z)$')
plt.ylabel(r'$\log$ SFRD $(M_\odot\ {\rm yr^{-1}\ Mpc^{-3}})$')
plt.legend()
plt.axis([-0.1,0.9,-2.15,-0.85])
plt.savefig('SFRDz.png',dpi=300)


plt.figure()
l2, = plt.plot(np.log10(1+z),np.log10(rho_g),'bo')
plt.plot(np.log10(1+z),np.log10(grho_g),'bd')

legend2 = plt.legend([l2],['gas'],loc=3)
ax = plt.gca().add_artist(legend2)
plt.xlabel(r'$\log\ (1+z)$')
plt.ylabel(r'$\log$ $\rho_{\rm met}$ $(M_\odot\ Mpc^{-3}})$')
plt.legend(loc = 4)
plt.axis([0,1,3,7])
plt.savefig('rhoz.png',dpi=300)

np.savez(p+p[:-1]+'_cosmo.npz',SFRD = SFRD, redshift = z,rho_g = rho_g, rho_g_gal = grho_g)
