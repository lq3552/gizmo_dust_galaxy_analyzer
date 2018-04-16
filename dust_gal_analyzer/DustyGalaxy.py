import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# TO DO: gas_Z be molecular gas_Z
class DustyGalaxy(object):

	def __init__(self, glist):
		self._Mg = glist['gas_mass'] # Msun
		self._Ms = glist['star_mass']
		self._Md = glist['dust_mass']
		self._Z  = glist['gas_Z']
		self._Zd = self._Md / self._Mg
		self._SFR =glist['SFR'] # Msun/yr
		self._SFRD = glist['SFRD'] # comoving SFR density, Msun/yr/Mpc^3
		self._z = glist['redshift']
		self._hp = glist['hp'] # Hubble Parameter
		self._dim = glist['dimension'] # comoving Mpc
#		self._dim = glist['dimension'] * (1+self._z) # comoving Mpc
		self._vol = (self._dim[0]*self._dim[1]*self._dim[2]) # comoving Mpc^3
		self._Od = np.sum(self._Md)/self._vol # only for test
		self._Og = np.sum(self._Z * self._Mg)/self._vol # only for test

		self._set_mpl()


	def get(self, field):
		if(field == 'gas_mass'): return self._Mg
		if(field == 'star_mass'): return self._Ms
		if(field == 'dust_mass'): return self._Md
		if(field == 'gas_Z'): return self._Z
		if(field == 'dust_Z'): return self._Zd
		if(field == 'SFR'): return self._SFR
		if(field == 'SFRD'): return self._SFRD
		if(field == 'redshift'): return self._z
		if(field == 'hubble_constant'): return self._hp
		if(field == 'dimension'): return self._dim
		if(field == 'volume'): 
			return self._vol
		else:
			return ['gas_mass','star_mass','dust_mass','gas_Z','dust_Z','SFR','SFRD','redshift','hubble_constant','dimension','volume']

########################################################################
#	routines for diagnostic plots
#
#	Note
#	---------------------------
#	All routines work like plt.plot plus plt.xlabel, plt.ylabel,
#	if you need to customize other elements 
#	of plots or save figures, you have to call 
#	other pyplot routines.
#
#	For example:
#	fig = plt.figure()
#	DustyGalaxy.plot_dmf(sty='r-',xlabel='M',ylabel=r'$\phi$')
#	plt.axis([0,8,-4.4,-1.0])
#	plt.savefig('DMF.png',dpi=300)
#########################################################################
	def plot_dmf(self,sty = 'o',label = None,yMIN = -1,yMAX = -1,\
			xlabel =r'$\log\ M_d\ [M_\odot]$',\
			ylabel =r'$\log\ \phi\ [{\rm Mpc^{-3}dex^{-1}}]$'):
		# plot dust mass function
		l1 = self._plot_dmf(sty,label,yMIN,yMAX,xlabel,ylabel)
		return l1
	
	def plot_gmf(self,sty = 'o',label = None,yMIN = -1,yMAX = -1,\
			xlabel =r'$\log\ M_g\ [M_\odot]$',\
			ylabel =r'$\log\ \phi\ [{\rm Mpc^{-3}dex^{-1}}]$'):
		# plot dust mass function
		l1 = self._plot_gmf(sty,label,yMIN,yMAX,xlabel,ylabel)
		return l1
	
	def plot_smf(self,sty = 'o',label = None,yMIN = -1,yMAX = -1,\
			xlabel =r'$\log\ M_*\ [M_\odot]$',\
			ylabel =r'$\log\ \phi\ [{\rm Mpc^{-3}dex^{-1}}]$'):
		# plot dust mass function
		l1 = self._plot_smf(sty,label,yMIN,yMAX,xlabel,ylabel)
		return l1
	
	def plot_dgr(self,x = None,sty = 'o',label = None,\
			xlabel = r'$12+\log(O/H)$',\
			ylabel = r'$\log$ (Gas/Dust mass ratio)' ):
		# plot dust to gas ratio (DGR) vs x (default: Mstar)
		if x == None:
			x = np.log10(self._Z/0.0134) + 8.69
		l1 = self._plot_dgr(x,sty,label,xlabel,ylabel)
		return l1

	def plot_dsr(self,x = None,sty = 'o',label = None,\
			xlabel = r'$\log\ M_*\ [M_\odot]$',\
			ylabel = r'$\log$ DSR' ):
		# plot dust to stellar mass ratio (DGR) vs x (default: Mstar)
		if x == None:
			x = self._Ms
		l1 = self._plot_dsr(x,sty,label,xlabel,ylabel)
		return l1

	def plot_dmr(self,x = None,sty = 'o',label = None,\
			xlabel = r'$\log\ Z\ [Z_\odot]$',\
			ylabel = r'$\log$ DMR' ):
		# plot dust to gas metallicity ratio (DMR) vs x (default: dust_Z)
		if x == None:
			x = self._Z/0.0134
		l1 = self._plot_dmr(x,sty,label,xlabel,ylabel)
		return l1

	def plot_mmr(self,sty = 'o',label = None,\
			xlabel = r'$\log\ (M_*/M_\odot)$',\
			ylabel= r'$\log\ (Z_g/Z_\odot)$'):
		l1 = self._plot_mmr(self._Ms,self._Z/0.0134,sty,label,xlabel,ylabel)
		return l1

	def _plot_dmf(self,sty,label,yMIN,yMAX,xlabel,ylabel):
		Md = self._Md[np.where(self._Md>0.0)].flatten()
		bins,mf = self._massfunc(Md,binsize = 0.1,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _plot_gmf(self,sty,label,yMIN,yMAX,xlabel,ylabel):
		Mg = self._Mg[np.where(self._Mg>0.0)].flatten()
		bins,mf = self._massfunc(Mg,binsize = 0.1,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _plot_smf(self,sty,label,yMIN,yMAX,xlabel,ylabel):
		Ms = self._Ms[np.where(self._Ms>0.0)].flatten()
		bins,mf = self._massfunc(Ms,binsize = 0.1,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _massfunc(self,Mass,binsize=0.2,MIN=-1,MAX=-1):
		"""Massfunction calculation

			Parameters
			----------
			binsize : float (optional)
			MIN : float (optional)
			MAX : float (optional)

			Returns
			-------
			bins : numpy_array
			massfnc : numpy_array
		"""
		vol = self._vol
		if MIN == -1:
			MIN = np.log10(np.amin(Mass))
		if MAX == -1:
			MAX = np.log10(np.amax(Mass))
		bins   = np.arange(MIN-1,MAX+1,binsize)
		elhist = np.histogram(np.log10(Mass),bins=bins)
		mf     = np.zeros(len(bins))
		for i in range(0,len(bins)-1):
			mf[i] = np.log10(float(elhist[0][i])/float(vol)/float(binsize)+1e-30)
		mf[len(bins)-1] = -30
		return (bins+0.5*binsize),mf
	
	def _plot_dgr(self,x,sty,label,xlabel,ylabel):
		l1, = plt.plot(x,-np.log10(self._Md/self._Mg),sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _plot_dsr(self,x,sty,label,xlabel,ylabel):
		l1, = plt.plot(np.log10(x),np.log10(self._Md/self._Ms),sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _plot_dmr(self,x,sty,label,xlabel,ylabel):
		l1, = plt.plot(np.log10(x),np.log10(self._Zd/self._Z),sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1
	
	def _plot_mmr(self,mass,met,sty,label,xlabel,ylabel):
		l1, = plt.plot(np.log10(mass),np.log10(met),sty,label=label)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _set_mpl(self):
		mpl.rcParams['figure.figsize'] = (10.0,8.0)
		mpl.rcParams['lines.linewidth'] = 3.0
		mpl.rcParams['axes.linewidth'] = 2.0
		mpl.rcParams['axes.labelsize'] = 24
		mpl.rcParams['xtick.labelsize'] = 22
		mpl.rcParams['ytick.labelsize'] = 22
		mpl.rcParams['xtick.major.width'] = 2.0
		mpl.rcParams['xtick.minor.width'] = 2.0
		mpl.rcParams['ytick.major.width'] = 2.0
		mpl.rcParams['ytick.minor.width'] = 2.0
		mpl.rcParams['legend.fontsize'] = 22
		mpl.rc('font',size=24)


if __name__ == "__main__":
	# usage example
	import sys
	if len(sys.argv) < 2: raise NameError('Please give the name of .npz file containing information of galaxies!')
	if sys.argv[1][-4:] != '.npz': raise NameError('The file name must end with .npz')

	data = np.load(sys.argv[1])
	gal =  DustyGalaxy(glist=data)

	print gal.get('')

	plt.figure()
	gal.plot_dmf()
	plt.axis([4,10,-4.4,-1.0])
	plt.savefig('DMF_'+sys.argv[1][:-4]+'.png',dpi=300)

	plt.figure()
	gal.plot_gmf()
	plt.axis([7,10.5,-4.4,-1.0])
	plt.savefig('GMF_'+sys.argv[1][:-4]+'.png',dpi=300)

	plt.figure()
	gal.plot_smf()
	plt.axis([8.5,11,-4.4,-1.0])
	plt.savefig('SMF_'+sys.argv[1][:-4]+'.png',dpi=300)

	plt.figure()
	gal.plot_dgr()
	plt.savefig('DGR_'+sys.argv[1][:-4]+'.png',dpi=300)

	plt.figure()
	gal.plot_dmr()
	plt.savefig('DMR_'+sys.argv[1][:-4]+'.png',dpi=300)

	plt.figure()
	gal.plot_dsr()
	plt.savefig('DSR_'+sys.argv[1][:-4]+'.png',dpi=300)
