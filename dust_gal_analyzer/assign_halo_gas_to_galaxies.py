from __future__ import print_function
from __future__ import division

import numpy as np

def periodic(x,halfbox,boxsize):
	if x < -halfbox:
		x += boxsize
	if x > halfbox:
		x -= boxsize
	return x

def assign_halo_gas_to_galaxies(
		# internal to the halo
		internal_galaxy_pos,
		internal_galaxy_mass,
		internal_glist,
		internal_galaxy_index_list,
		# global values
		Mg,
		galaxy_glist,
		gpos, gmass,
		boxsize, halfbox):
	ngas_internal = internal_glist.shape[0]
	n_internal_galaxies = internal_galaxy_mass.shape[0]

	for j in range(0,ngas_internal):
		i = internal_glist[j]
		gas_galaxy_index = galaxy_glist[i]
		max_index = 0
		max_mwd = 0.0

		if gas_galaxy_index < 0:
			for k in range(0,n_internal_galaxies):
				d2 = (periodic(gpos[i,0] - internal_galaxy_pos[k,0], halfbox, boxsize)**2 +
					  periodic(gpos[i,1] - internal_galaxy_pos[k,1], halfbox, boxsize)**2 +
					  periodic(gpos[i,2] - internal_galaxy_pos[k,2], halfbox, boxsize)**2)
				mwd = internal_galaxy_mass[k] / d2
				if mwd > max_mwd:
					max_mwd = mwd
					max_index = k

			gas_galaxy_index = internal_galaxy_index_list[max_index]
			Mg[gas_galaxy_index] += gmass[i]

	return Mg


def add_halo_gas_mass(obj,Mg,**kwargs):
# add halo gas mass (usually low density) to the galaxy gas mass which only considers SF gas
	boxsize = obj.simulation.boxsize.d
	halfbox = boxsize / 2.0
	galaxy_glist = np.array(obj.global_particle_lists.galaxy_glist,dtype=np.int32)

	# global gas properties
	gpos  = obj.data_manager.pos[obj.data_manager.glist]
	gmass = obj.data_manager.mass[obj.data_manager.glist]

	halos = obj.halos
	galaxy_pos = np.array([gal.pos for gal in obj.galaxies])
	galaxy_mass = np.array([gal.masses['total'] for gal in obj.galaxies])
	
	for h in halos:
		if len(h.galaxy_index_list) == 0:
			continue

		internal_galaxy_pos = galaxy_pos[h.galaxy_index_list]
		internal_galaxy_mass = galaxy_mass[h.galaxy_index_list]
		internal_galaxy_index_list = np.array(h.galaxy_index_list,dtype=np.int32)
		internal_glist = np.array(h.glist,dtype=np.int32)

		Mg = assign_halo_gas_to_galaxies(internal_galaxy_pos,
									internal_galaxy_mass,
									internal_glist,
									internal_galaxy_index_list,
									Mg,
									galaxy_glist,
									gpos, gmass,
									boxsize, halfbox)
	
	return Mg
