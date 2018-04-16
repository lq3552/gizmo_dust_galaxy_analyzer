import numpy as np
import os.path as path
import caesar as cs
import yt
import constants as C

class DustyGalaxyExtractor(object):
	###Extract galaxies and their properties from a GIZMO-dust snapshot
	###Parameter: 
	###    fname: name of snapshot
	###  replace: 0 -  if caesar output exists, directly load it and prohibit saving; 1 - not check if caesar output exists
	def __init__(self, fname, replace = 0):
		# load the snapshot into yt
		self.ds = yt.load(fname)
		self._file = fname
		self._replace = replace
		# load or calculate caesar object
		if (not (path.isfile('caesar_'+fname))) or (replace==1):
			# create a new CAESAR object, and pass along the yt dataset
			self.obj = cs.CAESAR(self.ds)
			# find haloes and galaxies
			self.obj.member_search()
		else:
			self.obj = cs.load('caesar_'+fname)
			self.obj.yt_dataset = self.ds

	def savec(self):
		# save the caesar output
		if (not (path.isfile('caesar_'+self._file))) or (self._replace==1):
			self.obj.save('caesar_'+self._file)
		else:
			print 'Saving aborted: caesar_'+self._file+' already exists!'
	
	def gal_extract(self):
		ds = self.ds
		obj = self.obj

		# properties of individual galaxies
		Mg = []
		Md = []
		Ms = []
		Z = []
		SFR = []
		for gal in obj.galaxies:
			index = gal.glist
			ad = ds.all_data()
			Mg.append(gal.masses['gas'].in_units('msun'))
			Ms.append(gal.masses['stellar'].in_units('msun'))
			Z.append(gal.metallicities['mass_weighted'])
			SFR.append(gal.sfr)
			Md.append(np.sum(ad[('PartType0', 'Dust_Masses')][index])*C.Mcode/C.Msun/ds.hubble_constant)
		Md = np.array(Md)
		Mg = np.array(Mg)
		Ms = np.array(Ms)
		Z = np.array(Z)
		SFR = np.array(SFR)
		
		# properties of the cosmo box
		dims = np.array(ds.domain_width.in_units('Mpccm'))
		SFRD = np.sum(ad[('PartType0', 'StarFormationRate')])/(dims[0]*dims[1]*dims[2])

		fname = 'gal_'+self._file.split('.')[0]+'.npz'
		np.savez(fname,gas_mass = Mg, dust_mass = Md, star_mass = Ms, gas_Z = Z, SFR = SFR,\
				hp = ds.hubble_constant, redshift = ds.current_redshift,dimension = dims,\
				SFRD = SFRD)
		print 'Save the table of galaxy properties to '+fname+'.'
		glist = np.load(fname)
		return glist



if __name__ == "__main__":
	# usage example
	import sys
	fname = sys.argv[1]
	dge = DustyGalaxyExtractor(fname,replace=0)
	print dge._file
	dge.savec()
	glist = dge.gal_extract()
	print glist['dust_mass'].shape,glist['star_mass'].shape
