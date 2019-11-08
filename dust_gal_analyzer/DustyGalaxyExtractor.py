from __future__ import division
from __future__ import print_function

import numpy as np
import os.path as path
import caesar as cs
import yt

class DustyGalaxyExtractor(object):
	###Extract galaxies and their properties from a GIZMO-dust snapshot
	###Parameter: 
	###    fname: name of snapshot
	###  replace: 0 -  if caesar output exists, directly load it and prohibit saving; 1 - not check if caesar output exists
	def __init__(self, snapshot, caesar, replace = 0, **kwargs):
		# paths, filenames and load snapshot
		self._yt = snapshot
		self._caesar = caesar
		self._replace = replace
		if (not (path.isfile(self._yt))): raise NameError('Snapshot file {} does not exist!'.format(self._yt))
		self.ds = yt.load(snapshot)
		_ = self._yt.rsplit('/',1)
		if len(_) == 1:
			self._dir = './'
		else:
			self._dir = _[0] + '/'

		if ('dust',True) in kwargs.items():
			self._IsDustActive = True
		else:
			self._IsDustActive = False

		# load or calculate caesar object
		if (not (path.isfile(self._caesar))) or (replace==1):
			# create a new CAESAR object, and pass along the yt dataset
			print('Create new CAESAR file...')
			self.obj = cs.CAESAR(self.ds)
			self.obj.member_search(**kwargs)#(blackholes=True,lowres=[2,3])
		else:
			self.obj = cs.load(self._caesar)
			self.obj.yt_dataset = self.ds

	def savec(self):
		# save the caesar output
		if (not (path.isfile(self._caesar))) or (self._replace==1):
			self.obj.save(self._caesar)
		else:
			print('Saving aborted: CAESAR file already exists!')
	
	def gal_extract(self):
		ds = self.ds
		obj = self.obj
		ad = ds.all_data()
		Mcode = ds.mass_unit.in_units('msun').value

		# properties of individual galaxies
		Mg,Mg_half,Md,Md_half,Ms,Mh,Mhi,Mh2,Mtot = [],[],[],[],[],[],[],[],[]
		rb,rbh,rg,rgh,rs,rsh,r200c,r200,rvir = [],[],[],[],[],[],[],[],[]
		Tmw,Tsw,Tvir = [],[],[]
		Z = []
		Zm = []
		Zs = []
		SFR = []
		for gal in obj.galaxies:
			rg2 = np.array([gal.radii['gas'].in_units('kpc')])[0]**2
			pos = np.array(gal.pos.in_units('kpc'))
			index = gal.glist
			indexs = gal.slist
			if (self._IsDustActive):
				indexd = gal.dlist
			data_Mg = np.array(ad[("PartType0","Masses")][index].in_units('msun'))
			if (self._IsDustActive):
				data_Md = np.array(ad[("PartType3","Masses")][indexd])
			else:
				data_Md = np.array(ad[("PartType0","Dust_Masses")][index])
			loc =  np.array(ad[("PartType0","Coordinates")][index].in_units('kpc'))
			loc = loc - pos
			rloc = loc[:,0]**2+loc[:,1]**2+loc[:,2]**2
			if (self._IsDustActive):
				locd = np.array(ad[("PartType3","Coordinates")][indexd].in_units('kpc'))
				locd = locd - pos
				rlocd = locd[:,0]**2+locd[:,1]**2+locd[:,2]**2
			filt = np.where(rloc <= 0.25*rg2)
			if (self._IsDustActive):
				filtd = np.where(rlocd <= 0.25*rg2)
			else:
				filtd = filt
			Mg.append(gal.masses['gas'].in_units('msun'))
			Ms.append(gal.masses['stellar'].in_units('msun'))
			Mh.append(gal.masses['H'].in_units('msun'))
			Mhi.append(gal.masses['HI'].in_units('msun'))
			Mh2.append(gal.masses['H2'].in_units('msun'))
			Mtot.append(gal.masses['total'].in_units('msun'))
			if (self._IsDustActive):
				Md.append(np.sum(ad[('PartType3', 'Masses')][indexd])*Mcode)
			else:
				Md.append(np.sum(ad[('PartType0', 'Dust_Masses')][index])*Mcode)
			Mg_half.append(np.sum(data_Mg[filt]))
			Md_half.append(np.sum(data_Md[filtd])*Mcode)
			rb.append(gal.radii['baryon'].in_units('kpc'))
			rbh.append(gal.radii['baryon_half_mass'].in_units('kpc'))
			rg.append(gal.radii['gas'].in_units('kpc'))
			rgh.append(gal.radii['gas_half_mass'].in_units('kpc'))
			rs.append(gal.radii['stellar'].in_units('kpc'))
			rsh.append(gal.radii['stellar_half_mass'].in_units('kpc'))
			r200c.append(gal.radii['r200c'].in_units('kpc'))
			r200.append(gal.radii['r200'].in_units('kpc'))
			rvir.append(gal.radii['virial'].in_units('kpc'))
			Tvir.append(gal.temperatures['virial'].in_units('K'))
			Tmw.append(gal.temperatures['mass_weighted'].in_units('K'))
			Tsw.append(gal.temperatures['sfr_weighted'].in_units('K'))
			Z.append(gal.metallicities['sfr_weighted'])
			Zm.append(gal.metallicities['mass_weighted'])
			w =  np.sum(np.array(ad[('PartType4', 'Masses')][indexs]).flatten())
			if len(np.array(ad[('PartType4', 'Masses')][indexs]).flatten()) > 0:
				if w > 0.:
					Zs.append(np.average(np.array(ad[('PartType4', 'Metallicity_00')][indexs]).flatten(),\
								weights = np.array(ad[('PartType4', 'Masses')][indexs]).flatten()))
				else:
					Zs.append(np.average(np.array(ad[('PartType4', 'Metallicity_00')][indexs]).flatten(),\
								weights = np.array(ad[('PartType4', 'Masses')][indexs]).flatten()))
			else:
				Zs.append(0.0)
			SFR.append(gal.sfr)
		Md = np.array(Md)
		Mg = np.array(Mg)
		Md_half = np.array(Md)
		Mg_half = np.array(Mg)
		if not (self._IsDustActive):
			Mg -= Md
			Mg_half -= Md_half
		Ms = np.array(Ms)
		Mh = np.array(Mh)
		Mhi = np.array(Mhi)
		Mh2 = np.array(Mh2)
		Mtot = np.array(Mtot)
		rb = np.array(rb)
		rbh = np.array(rbh)
		rg = np.array(rg)
		rgh = np.array(rgh)
		rs = np.array(rs)
		rsh = np.array(rsh)
		r200 = np.array(r200)
		r200c = np.array(r200c)
		rvir = np.array(rvir)
		Tmw  = np.array(Tmw)
		Tsw  = np.array(Tsw)
		Tvir  = np.array(Tvir)
		Z = np.array(Z)
		Zs = np.array(Zs)
		SFR = np.array(SFR)
		
		# properties of the cosmo box
		dims = np.array(ds.domain_width.in_units('Mpccm'))
		SFRD = np.sum(ad[('PartType0', 'StarFormationRate')])/(dims[0]*dims[1]*dims[2])
		if (self._IsDustActive):
			rhod = np.sum(ad[('PartType3', 'Masses')].in_units('msun'))/(dims[0]*dims[1]*dims[2])
		else:
			rhod = np.sum(ad[('PartType0', 'Dust_Masses')])*Mcode/(dims[0]*dims[1]*dims[2])
		massg = np.array(ad[('PartType0', 'Masses')].in_units('msun'))
		if not (self._IsDustActive):
			massg -= np.array(ad[('PartType0', 'Dust_Masses')] * Mcode)
		rhog = np.sum(massg)/(dims[0]*dims[1]*dims[2]) 
		rhogz = np.sum(massg*ad[('PartType0','Metallicity_00')])/(dims[0]*dims[1]*dims[2])

		fname = self._dir + 'gal_'+self._yt.rsplit('/',1)[-1].split('.')[0]+'.npz'
		np.savez(fname,\
				gas_mass = Mg, dust_mass = Md,\
				gas_mass_hr = Mg_half, dust_mass_hr = Md_half,\
				star_mass = Ms, star_Z = Zs, gas_Z = Z,gas_Zm = Zm, SFR = SFR,\
				mass_h = Mh, mass_hi = Mhi, mass_h2 = Mh2, mass_tot = Mtot,\
				baryon_radii = rb,baryon_radii_hm = rbh,\
				star_radii = rs,star_radii_hm = rsh,\
				gas_radii = rg,gas_radii_hm = rgh,\
				r200 = r200,r200c = r200c, radii_vir = rvir,\
				T_mass_weighted = Tmw, T_sfr_weighted = Tsw, T_vir = Tvir,\
				hp = ds.hubble_constant, redshift = ds.current_redshift,dimension = dims,\
				SFRD = SFRD,rhod = rhod,rhog = rhog,rhogz = rhogz)
		print('Save the table of galaxy properties to '+fname+'.')
		glist = np.load(fname)
		return glist



if __name__ == "__main__":
	# usage example
	import sys
	fname = sys.argv[1]
	dge = DustyGalaxyExtractor(fname,replace=0,blackholes=True,lowres=[2,3])
	print(dge._file)
	dge.savec()
	glist = dge.gal_extract()
