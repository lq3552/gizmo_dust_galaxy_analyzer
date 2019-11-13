from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

## TODO: use kwargs for pyplot paramaters

class DustyGalaxy(object):

	def __init__(self, glist):

		if glist[-4:] != '.npz': raise NameError('The filename must have extension .npz')

		self._glist = glist
		glist = np.load(glist)
		self._Mg = glist['gas_mass'] # Msun
		self._Mghi = glist['mass_hi']
		self._Mgh2 = glist['mass_h2']
		self._Mghr = glist['gas_mass_hr']
		self._Ms = glist['star_mass']
		self._Md = glist['dust_mass']
		self._Mdhr = glist['dust_mass_hr']
		self._rg = glist['gas_radii']
		self._Zg = glist['gas_Z']
		self._Zgm = glist['gas_Zm']
		if 'star_Z' in glist.files: # to do: a smarter way to gal.npz version backward compatibility
			self._Zs = glist['star_Z']
		else:
			self._Zs = 0.
		self._Zd = self._Md / glist['gas_mass']
		self._Ztot = (self._Zg * self._Mg + self._Md)/glist['gas_mass'] # depreciated?
		self._SFR =glist['SFR'] # Msun/yr
		self._SFRD = glist['SFRD'] # comoving SFR density, Msun/yr/Mpc^3
		self._rhod = glist['rhod'] # comoving dust density, Msun/Mpc^3
		self._rhog = glist['rhog']
		self._rhogz = glist['rhogz']
		self._z = glist['redshift']
		self._hp = glist['hp'] # Hubble Parameter
		self._dim = glist['dimension'] # comoving Mpc
		self._vol = (self._dim[0]*self._dim[1]*self._dim[2]) # comoving Mpc^3

		self._set_mpl()

	def get(self, field): # cumbersome, need improvement
		if(field == 'gas_mass'): return self._Mg
		if(field == 'gas_mass_hr'): return self._Mghr
		if(field == 'gas_radii'): return self._rg
		if(field == 'gas_Z'): return self._Zg
		if(field == 'gas_Zm'): return self._Zgm
		if(field == 'star_mass'): return self._Ms
		if(field == 'star_Z'): return self._Zs
		if(field == 'dust_mass'): return self._Md
		if(field == 'dust_mass_hr'): return self._Mdhr
		if(field == 'dust_Z'): return self._Zd
		if(field == 'dust_density'): return self._rhod
		if(field == 'SFR'): return self._SFR
		if(field == 'SFRD'): return self._SFRD
		if(field == 'gas_density'): return self._rhog
		if(field == 'gas_metal_density'): return self._rhogz
		if(field == 'redshift'): return self._z
		if(field == 'hubble_constant'): return self._hp
		if(field == 'dimension'): return self._dim
		if(field == 'volume'): 
			return self._vol
		else:
			return ['gas_mass','gas_mass_hr','gas_radii','gas_Z','gas_Zm',
			'star_mass','star_Z',
			'dust_mass','dust_mass_hr','dust_Z',
			'SFR',
			'SFRD','dust_density','gas_density','gas_metal_density','redshift','hubble_constant',
			'dimension','volume']

########################################################################
#	Analysis tools
########################################################################
	def get_ms_and_quenched_galaxies(self,eps=0.1,min_samples=30):
	# For now it generates a new npz and reload it to return a new DustyGalaxy object. I need to find a way to do this without S/L external files
		from sklearn.cluster import DBSCAN

		ssfr = self._SFR/self._Ms*1.e9
		Zmet = np.log10(self._Zg/0.0134)
		dgr = np.log10(self._Md/self._Mg)
		ssfr = np.log10(ssfr+10**(-2.5+0.3*self._z))
		index  = np.arange(len(Zmet))

		X = np.zeros((len(Zmet),4))
		X[:,0] = Zmet
		X[:,1] = ssfr 
		X[:,2] = dgr
		X[:,3] = index

		# clean up data and prepare vectors needed for fitting
		X = X[~np.isnan(X).any(axis=1)]
		X = X[~np.isinf(X).any(axis=1)]
		X = X[np.lexsort(np.rot90(X))] # sort the data in lexicographical order
		Zmet = X[:,0]
		dgr = X[:,2]
		index = X[:,3]

		# ND points for clustering (seperate "main sequence" from quenched galaxies)
		NN = len(X[:,0])
		XC = np.zeros((NN,3))
		XC[:,0] = X[:,0]
		XC[:,1] = X[:,1]
		XC[:,2] = X[:,2]

		# clustering
		clt = DBSCAN(eps=eps,min_samples=min_samples)
		clt.fit(XC)

		index_ms       = index[np.where(clt.labels_==0)].astype(int)
		index_quenched = index[np.where(clt.labels_!=0)].astype(int)
		glist_ms       = self._save_glist_sub(index_ms,fname='tmp_sf.npz')
		glist_quenched = self._save_glist_sub(index_quenched,fname='tmp_q.npz')

		return DustyGalaxy(glist_ms),DustyGalaxy(glist_quenched) 

	def _save_glist_sub(self,index,fname='tmp.npz'):
		glist = self._glist
		np.savez(fname,gas_mass=glist['gas_mass'][index],
			gas_mass_hr = glist['gas_mass_hr'][index],
			mass_hi = glist['mass_hi'][index],
			mass_h2 = glist['mass_h2'][index],
			star_mass = glist['star_mass'][index],
			dust_mass = glist['dust_mass'][index],
			dust_mass_hr = glist['dust_mass_hr'][index],
			gas_radii = glist['gas_radii'][index],
			gas_Z = glist['gas_Z'][index],
			SFR =glist['SFR'][index],
			SFRD = glist['SFRD'],
			rhod = glist['rhod'],
			rhog = glist['rhog'],
			rhogz = glist['rhogz'],
			redshift = glist['redshift'],
			hp = glist['hp'],
			dimension = glist['dimension'])
		return np.load(fname)
		


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
	def plot_dmf(self,binsize = 0.1,yMIN = -1,yMAX = -1,\
			sty = 'o',\
			xlabel =r'$\log\ M_d\ [M_\odot]$',\
			ylabel =r'$\log\ \phi\ [{\rm Mpc^{-3}dex^{-1}}]$',**kwargs):
		# plot dust mass function
		l1 = self._generic_plot_massfunc(self._Md,binsize,yMIN,yMAX,sty,xlabel,ylabel,**kwargs)
		return l1
	
	def plot_gmf(self,binsize = 0.1,yMIN = -1,yMAX = -1,\
			sty = 'o',\
			xlabel =r'$\log\ M_g\ [M_\odot]$',\
			ylabel =r'$\log\ \phi\ [{\rm Mpc^{-3}dex^{-1}}]$',**kwargs):
		# plot dust mass function
		l1 = self._generic_plot_massfunc(self._Mg,binsize,yMIN,yMAX,sty,xlabel,ylabel,**kwargs)
		return l1
	
	def plot_smf(self,binsize = 0.1,yMIN = -1,yMAX = -1,\
			sty = 'o',\
			xlabel =r'$\log\ M_*\ [M_\odot]$',\
			ylabel =r'$\log\ \phi\ [{\rm Mpc^{-3}dex^{-1}}]$',**kwargs):
		# plot stellar mass function
		l1 = self._generic_plot_massfunc(self._Ms,binsize,yMIN,yMAX,sty,xlabel,ylabel,**kwargs)
		return l1
	
	def plot_dgr(self,x = None, mode = 'scatter',\
			xlabel = r'$\log\ (Z/Z_\odot)$',\
			ylabel = r'$\log$ DGR', **kwargs):
		# plot dust to gas ratio (DGR) vs x (default: Zg)
		if x == None:
			x = np.log10(self._Zg/0.0134)
		l1 = self._generic_plot_scaling(x,np.log10(self._Md/self._Mg),mode,xlabel,ylabel,**kwargs)
		return l1

	def plot_dtm(self,x = None, mode = 'scatter',\
			xlabel = r'$\log\ (Z/Z_\odot)$',\
			ylabel = r'$\log$ DTM', **kwargs):
		# plot dust to (gas) metal ratio (DTM) vs x (default: Zg)
		if x == None:
			x = np.log10(self._Zg/0.0134)
		l1 = self._generic_plot_scaling(x,np.log10(self._Md/(self._Zg*self._Mg)),mode,xlabel,ylabel,**kwargs)
		return l1

	def plot_dsr(self,x = None, mode = 'scatter',\
			xlabel = r'$\log\ M_*\ [M_\odot]$',\
			ylabel = r'$\log$ DSR', **kwargs):
		# plot dust to stellar mass ratio (DGR) vs x (default: Mstar)
		if x == None:
			x = np.log10(self._Ms)
		l1 = self._generic_plot_scaling(x,np.log10(self._Md/self._Ms),mode,xlabel,ylabel,**kwargs)
		return l1

	def plot_mzr(self, y = None, mode = 'scatter',\
			xlabel = r'$\log\ (M_*/M_\odot)$',\
			ylabel= r'$\log\ Z\ [Z_\odot]$', **kwargs):
		# plot mass-metallicity relations
#		l1 = self._plot_mmr(self._Ms,np.log10(self._Zg/0.0134) + 8.69,sty,label,alpha,xlabel,ylabel)
#		l1 = self._plot_mmr(self._Ms,np.log10(self._ZgO/8.65e-3) + 8.69 + 0.4,sty,label,alpha,xlabel,ylabel)
		if y == None:
			y = np.log10(self._Zg/0.0134)
		l1 = self._generic_plot_scaling(np.log10(self._Ms),y,mode,xlabel,ylabel,**kwargs)
		return l1

	def plot_sdr(self,sty = 'o',label = None,\
			alpha = None,\
			xlabel = r'$\log\ (M_*/M_\odot)$',\
			ylabel= r'$\log\ (M_d/M_\odot)$'):
		l1 = self._plot_sdr(self._Ms,self._Md,sty,label,alpha,xlabel,ylabel)
		return l1

	def _plot_dmf(self,yMIN,yMAX,sty,xlabel,ylabel,**kwargs):
		Md = self._Md[np.where(self._Md>0.0)].flatten()
		bins,mf = self._massfunc(Md,binsize = 0.1,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,**kwargs)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _plot_gmf(self,yMIN,yMAX,sty,xlabel,ylabel,**kwargs):
		Mg = self._Mg[np.where(self._Mg>0.0)].flatten()
		bins,mf = self._massfunc(Mg,binsize = 0.1,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,**kwargs)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _plot_smf(self,yMIN,yMAX,sty,xlabel,ylabel,**kwargs):
		Ms = self._Ms[np.where(self._Ms>0.0)].flatten()
		bins,mf = self._massfunc(Ms,binsize = 0.1,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,alpha=alpha,label=label)
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
	
	def _generic_plot_massfunc(self,M,binsize,yMIN,yMAX,sty,xlabel,ylabel,**kwargs):
		M = M[np.where(M>0.0)].flatten()
		bins,mf = self._massfunc(M,binsize=binsize,MIN=yMIN,MAX=yMAX)
		l1, = plt.plot(bins,mf,sty,**kwargs)
		plt.xlabel(xlabel)
		plt.ylabel(ylabel)
		return l1

	def _generic_plot_scaling(self,x,y,mode,xlabel,ylabel,**kwargs):
		if mode == 'scatter':
			l1 = plt.scatter(x,y,**kwargs)
		elif mode == 'hexbin':
			l1 = plt.hexbin(x,y,**kwargs)
		else:
			raise AttributeError('Only "scatter" and "hexbin" are allowed!')
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
	glist = sys.argv[1]
	gal =  DustyGalaxy(glist)

	print(gal.get(''))

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
