import h5py
import os
import glob
import re
import numpy as np
from . import peano
import warnings
from scipy.integrate import quad

base_path = os.environ['EAGLE_BASE_PATH']
release = os.environ['EAGLE_ACCESS_TYPE']

class Snapshot:
    """ Basic SnapShot superclass which finds the relevant files and gets relevant information
    regarding the snapshot specified.

    arguments:
        run - the run (e.g. L0012N0188)
        model - an EAGLE model (e.g. Ref)
        tag - a tag string specifying a snapshot output (e.g. 028_z000p000)

    history:
        written - Mackereth (UoB) - 22/11/2019
     """
    def __init__(self, run, model, tag, load_particles=False):
        #store the snapshot identity info
        self.run = run
        self.model = model
        self.tag = tag
        if release == 'public':
            self.simlabel = self.model+self.run
            self.snaplabel = 'snapshot_'+self.tag
            self.base_subfile = 'snap_'+self.tag
            self.path = os.path.join(base_path, self.simlabel, self.snaplabel)
        elif release == 'ARI':
            self.snaplabel = 'snapshot_'+self.tag
            self.base_subfile = 'snap_'+self.tag
            self.path = os.path.join(base_path, self.run, self.model, 'data', self.snaplabel)
        else:
            raise Exception('private/custom data access is not yet implemented!')
        if not os.path.exists(os.path.join(self.path, self.base_subfile+'.0.hdf5')):
            raise Exception('could not see snapshot data in directory: '+self.path)
        #get the files related to this snapshot and load some of their metadata
        self.files = natural_sort(glob.glob(os.path.join(self.path, self.base_subfile+'*.hdf5')))
        self.nfiles = len(self.files)
        self.header_dict = dict(h5py.File(self.files[0], 'r')['/Header'].attrs.items())
        self.abundance_dict = dict(h5py.File(self.files[0], 'r')['/Parameters/ChemicalElements'].attrs.items())
        self.elements  = ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen', 'Silicon', 'Sulphur', 'Magnesium', 'Iron']
        self.solar_abundances = dict([(self.elements[i],self.abundance_dict['SolarAbundance_%s' % self.elements[i]]) for i in range(len(self.elements))])
        self.BoxSize = self.header_dict['BoxSize']
        self.HubbleParam  = self.header_dict['HubbleParam']
        self.Omega0, self.OmegaLambda, self.OmegaBaryon, self.a0 = self.header_dict['Omega0'], self.header_dict['OmegaLambda'], self.header_dict['OmegaBaryon'], self.header_dict['ExpansionFactor']
        self.NumPartTotal = self.header_dict['NumPart_Total']
        self.ParticleTypes = np.array([0,1,2,3,4,5])
        self.ParticleTypePresent = self.NumPartTotal > 0
        self.ParticleTypePresent_file = np.zeros((len(self.files),len(self.NumPartTotal)), dtype=bool)
        for ii, file in enumerate(self.files):
            head = dict(h5py.File(file, 'r')['/Header'].attrs.items())
            self.ParticleTypePresent_file[ii, head['NumPart_ThisFile'] > 0] = True
        self._ptypeind = {self.ParticleTypes[self.ParticleTypePresent][i]:i for i in range(len(self.ParticleTypes[self.ParticleTypePresent]))}
        #get the Hash Table info for P-H key sorting
        self.HashBits = dict(h5py.File(self.files[0], 'r')['/HashTable'].attrs.items())['HashBits']
        self.HashGridSideLength = 2**self.HashBits
        self.HashGridCellSize = self.BoxSize/self.HashGridSideLength
        self.firstkeys = np.zeros((len(self.ParticleTypes[self.ParticleTypePresent]),self.nfiles))
        self.lastkeys  = np.zeros((len(self.ParticleTypes[self.ParticleTypePresent]),self.nfiles))
        self.datasets = {}
        for ii,parttype in enumerate(self.ParticleTypes[self.ParticleTypePresent]):
            self.firstkeys[ii] = np.array(h5py.File(self.files[0], 'r')['/HashTable/PartType'+str(parttype)+'/FirstKeyInFile'])
            self.lastkeys[ii] = np.array(h5py.File(self.files[0], 'r')['/HashTable/PartType'+str(parttype)+'/LastKeyInFile'])
            #be sure we get a file with this parttype (only really an issue for when low N stars!!)
            ind = np.nonzero(h5py.File(self.files[0], 'r')['/HashTable/PartType'+str(parttype)+'/LastKeyInFile'][:])[0][0]
            self.datasets['PartType'+str(parttype)] = list(h5py.File(self.files[ind], 'r')['/PartType'+str(parttype)].keys())
        if load_particles:
            self._get_coordinates()

    def _get_coordinates(self):
        """ Load all the coordinates of the available particles
        """
        #load coordinates and velocities
        coordinates = []
        velocities = []
        for ii,type in enumerate(self.ParticleTypes[self.ParticleTypePresent]):
            #now load the coordinates in these files and save the indices for each particle type
            thistypecoord, thistypevels = self._get_parttype_indices(type, self.files)
            coordinates.append(thistypecoord)
            velocities.append(thistypevels)
        self.velocities = velocities
        self.coordinates = coordinates

    def _get_parttype_indices(self, parttype, files):
        """get the coordinates and indices for a given particle type in a given region"""
        coords, velocities, indices = [], [], []
        for ii,file in enumerate(files):
            #check this particle type is present here
            if not _particle_type_present(parttype, file):
                return None, None
            # load the file
            thisfilecoords = np.array(h5py.File(file, 'r')['/PartType'+str(parttype)+'/Coordinates'])
            thisfilevels = np.array(h5py.File(file, 'r')['/PartType'+str(parttype)+'/Velocity'])
            #store the coordinates and the indices of these particles in the file
            coords.append(thisfilecoords)
            velocities.append(thisfilevels)
        return np.concatenate(coords), np.concatenate(velocities)

    def _get_coords_vels(self, parttype, files):
        """get the coordinates and velocities for all particles of a certain type"""
        if not self.ParticleTypePresent[parttype]:
            warnings.warn('Particle type is not present, returning empty arrays...')
            return np.array([]), np.array([]), np.array([])
        coords, velocities, indices = [], [], []
        for file in files:
            # load the file
            thisfilecoords = np.array(h5py.File(file, 'r')['/PartType'+str(parttype)+'/Coordinates'])
            thisfilevels = np.array(h5py.File(file, 'r')['/PartType'+str(parttype)+'/Velocity'])
            #store the coordinates and the indices of these particles in the file
            coords.append(thisfilecoords)
            velocities.append(thisfilevels)
        return np.concatenate(coords), np.concatenate(velocities)

    def get_dataset(self, parttype, dataset, physical=False, cgs=False):
        """ get the data for a given entry in the HDF5 file for the given region """
        if not self.ParticleTypePresent[parttype]:
            warnings.warn('Particle type is not present, returning empty arrays...')
            return np.array([])
        key = os.path.join('/PartType'+str(parttype),dataset)
        if physical:
            #find conversion factor
            factor = self._conversion_factor(key, self.a0, self.HubbleParam, cgs=cgs)
        elif not physical and cgs:
            factor = h5py.File(self.files[0], 'r')[key].attrs['CGSConversionFactor']
        else:
            #else just multiply by 1!
            factor = 1.
        out = []
        for ii,file in enumerate(self.files):
            # load this file and get the particles
            out.append(np.array(h5py.File(file, 'r')[key]) * factor)
        return np.concatenate(out)

    def _conversion_factor(self, key, a, h, cgs=False):
        aexp_scale, h_scale = self._get_conversion_factor_exponents(key)
        if cgs:
            cgs_factor = h5py.File(self.files[0], 'r')[key].attrs['CGSConversionFactor']
        else:
            cgs_factor = 1.
        return a**(aexp_scale)*h**(h_scale)*cgs_factor

    def _get_conversion_factor_exponents(self, key):
        aexp_scale = h5py.File(self.files[0], 'r')[key].attrs['aexp-scale-exponent']
        h_scale = h5py.File(self.files[0], 'r')[key].attrs['h-scale-exponent']
        return aexp_scale, h_scale

    def _single_X_H(self,X,H,element):
        solar = self.solar_abundances[element]
        solarH = self.solar_abundances['Hydrogen']
        return np.log10(X/H)-np.log10(solar/solarH)

    def abundance_ratios(self,gas=False,smoothed=True):
        """ Compute element abundance ratios for the region, returns a dict of [X/H] """
        if smoothed:
            e_key = 'SmoothedElementAbundance'
        else:
            e_key = 'ElementAbundance'
        if gas:
            parttype = 0
        else:
            parttype = 4
        entries = []
        H = self.get_dataset(parttype,os.path.join(e_key,'Hydrogen'))
        for i in range(len(self.elements)):
            if self.elements[i] == 'Hydrogen' or self.elements[i] == 'Sulphur':
                continue
            X = self.get_dataset(parttype,os.path.join(e_key,self.elements[i]))
            entries.append((self.elements[i],self._single_X_H(X,H,self.elements[i])))
        return dict(entries)

    def t_lookback(self,a):
        return a / (np.sqrt(self.Omega0 * a + self.OmegaLambda * (a ** 4)))

    def z2age(self,z):
        a = 1 / (1 + z)
        t = np.array([quad(self.t_lookback, x, self.a0)[0] for x in a])
        return (1 / (self.HubbleParam * 100)) * (3.086e19 / 3.1536e16) * t

    def a2age(self,a):
        t = np.array([quad(self.t_lookback, x, self.a0)[0] for x in a])
        return (1 / (self.HubbleParam * 100)) * (3.086e19 / 3.1536e16) * t

    def z2tau(self,z):
        t_em = quad(self.t_lookback, 0., self.a0)[0]
        t_em = (1 / (self.HubbleParam * 100)) * (3.086e19 / 3.1536e16) * t_em
        a = 1 / (1 + z)
        t = np.array([quad(self.t_lookback, x, self.a0)[0] for x in a])
        return t_em - ((1 / (self.HubbleParam * 100)) * (3.086e19 / 3.1536e16) * t)

    def a2tau(self,a):
        t_em = quad(self.t_lookback, 0., self.a0)[0]
        t_em = (1 / (self.HubbleParam * 100)) * (3.086e19 / 3.1536e16) * t_em
        t = np.array([quad(self.t_lookback, x, self.a0)[0] for x in a])
        return t_em - ((1 / (self.HubbleParam * 100)) * (3.086e19 / 3.1536e16) * t)






class SnapshotRegion(Snapshot):
    """ A class inheriting from SnapShot, which defines a region inside a larger simulation snapshot.
    when initialised, this will read the files in that region, and get the indices of the particles inside the
    desired region. The necessary datasets can then be loaded by using get_dataset.

    arguments:
        run - the run (e.g. L0012N0188)
        model - an EAGLE model (e.g. Ref)
        tag - a tag string specifying a snapshot output (e.g. 028_z000p000)
        center - the center of the desired region
        sidelength - the length of a side of the volume required

    history:
        written - Mackereth (UoB) - 22/11/2019
     """
    def __init__(self, run, model, tag, center, sidelength, just_get_files=False):
        #we want everything from SnapShot plus some extras
        super().__init__(run, model, tag)
        self.center = center
        self.sidelength = sidelength
        self.centered = False
        self._index_region(self.center, self.sidelength, justfiles=just_get_files)

    def _index_region(self, center, side_length, phgrid_n=70, justfiles=False):
        """ Load a region defined by a central cordinate and a side length
        arguments:
        center - the [x,y,z] coordinate of the desired center (simulation units)
        side_length - the desired side length (in the simulation units)

        keyword arguments:
        phgrid_n - the number of grid points along a side length to look for PH cells (default 70)
        """
        #work out which files contain the desired region
        grid = peano.coordinate_grid(center, side_length, self.BoxSize, n=phgrid_n)
        keys  = peano.get_unique_grid_keys(grid, self.HashGridCellSize, self.BoxSize, bits=self.HashBits)
        particles_in_volume = self.ParticleTypes[self.ParticleTypePresent]
        self.files_for_region = []
        self.file_indices = []
        coordinates = []
        velocities = []
        indices = []
        for ii in self.ParticleTypes:
            if not self.ParticleTypePresent[ii]:
                continue
            Nfiles = self._get_parttype_files(ii, keys)
            if len(Nfiles) < 1:
                #particle is not present in the region - remove from here
                self.ParticleTypePresent[ii]  = 0
                continue
            thisfiles = np.array(self.files)[Nfiles]
            thisindices = Nfiles
            self.files_for_region.append(thisfiles)
            self.file_indices.append(Nfiles)
            if justfiles:
                continue
            present = False
            for file in thisfiles:
                present += _particle_type_present(ii, file)
            if present:
                #now load the coordinates in these files and save the indices for each particle type
                thistypecoord, thistypevels, thistypeindices = self._get_parttype_indices(ii, thisfiles, thisindices)
                if thistypecoord is None:
                    self.ParticleTypePresent[ii] = 0
                    continue
                coordinates.append(thistypecoord)
                velocities.append(thistypevels)
                indices.append(thistypeindices)
            else:
                self.ParticleTypePresent[ii] = 0
        if not justfiles:
            self.velocities = velocities
            self.coordinates = coordinates
            self.indices = indices
            self.NumPart_ThisRegion = np.zeros(len(self.NumPartTotal),dtype=np.int64)
            for ii,type in enumerate(self.ParticleTypes[self.ParticleTypePresent]):
                self.NumPart_ThisRegion[type] = len(self.coordinates[ii])


    def _get_parttype_indices(self, parttype, files, file_indices):
        """get the coordinates and indices for a given particle type in a given region"""
        coords, velocities, indices = [], [], []
        for ii,file in enumerate(files):
            #check this particle type is present here
            if not _particle_type_present(parttype, file):
                return None, None, None
            # load the file
            thisfilecoords = np.array(h5py.File(file, 'r')['/PartType'+str(parttype)+'/Coordinates'])
            thisfilevels = np.array(h5py.File(file, 'r')['/PartType'+str(parttype)+'/Velocity'])
            if (np.array(self.center)+self.sidelength > self.BoxSize).any():
                thisfilecoords = thisfilecoords - (self.center - self.BoxSize/2.)
                thisfilecoords = np.mod(thisfilecoords,self.BoxSize)
                thisfilecoords -= self.BoxSize/2.
                thisfilecoords += self.center
            # mask it to the region desired
            mask = (np.fabs(thisfilecoords[:,0]-self.center[0]) < self.sidelength/2.) &\
                    (np.fabs(thisfilecoords[:,1]-self.center[1]) < self.sidelength/2.) &\
                     (np.fabs(thisfilecoords[:,2]-self.center[2]) < self.sidelength/2.)
            #store the coordinates and the indices of these particles in the file
            thisfileindices = np.where(mask)[0]
            coords.append(thisfilecoords[mask])
            velocities.append(thisfilevels[mask])
            indices.append(thisfileindices)
        return np.concatenate(coords), np.concatenate(velocities), indices

    def get_dataset(self, parttype, dataset, physical=False, cgs=False):
        """ get the data for a given entry in the HDF5 file for the given region """
        if not self.ParticleTypePresent[parttype]:
            warnings.warn('Particle type is not present, returning empty arrays...')
            return np.array([])
        key = os.path.join('/PartType'+str(parttype),dataset)
        if physical:
            #find conversion factor
            factor = self._conversion_factor(key, self.a0, self.HubbleParam, cgs=cgs)
        elif not physical and cgs:
            factor = h5py.File(self.files[0], 'r')[key].attrs['CGSConversionFactor']
        else:
            #else just multiply by 1!
            factor = 1.
        out = []
        ptypeind = self._ptypeind[parttype]
        for ii,file in enumerate(self.files_for_region[ptypeind]):
            if not _particle_type_present(parttype, file):
                continue
            # load this file and get the particles
            out.append(np.array(h5py.File(file, 'r')[key])[self.indices[ptypeind][ii]] * factor)
        if len(out) < 2:
            return out[0]
        return np.concatenate(out)

    def _get_parttype_files(self, parttype, keys):
        """ get the files containing this region for a given particle type """
        Nfiles = []
        ptypeind = self._ptypeind[parttype]
        for i in range(len(keys)):
            if len(np.where(self.firstkeys[ptypeind] < keys[i])[0]) < 1:
                start = 0
            else:
                start = np.where(self.firstkeys[ptypeind] < keys[i])[0][-1]
            if len(np.where(self.firstkeys[ptypeind] > keys[i])[0]) < 1:
                end = len(self.firstkeys[ptypeind])-1
            else:
                end = np.where(self.firstkeys[ptypeind] > keys[i])[0][0]
            Nfiles.extend(np.arange(start,end+1,1))
        Nfiles = np.unique(Nfiles)
        return Nfiles


    def angular_momentum(self, parttype=1, percentile=10):
        """Compute the angular momentum of particles within some percentile of their
           radii
        """
        ptypeind = self._ptypeind[parttype]
        pos, vel = self.coordinates[ptypeind], self.velocities[ptypeind]
        radii = np.linalg.norm(self.coordinates[ptypeind], axis=1)
        inside = radii < np.percentile(radii, percentile)
        if parttype == 1:
            mass= np.ones(len(pos))*self.header_dict['MassTable'][1]
        else:
            mass = self.get_dataset(parttype, 'Mass')
        vec = np.cross(pos[inside],vel[inside]*mass[inside][:,np.newaxis])
        tot = np.sum(vec, axis=0)
        return tot/np.linalg.norm(tot)

    def _transform(self,vector):
        """Build a transformation matrix"""
        a = vector
        b = np.matrix([0,0,1])
        v = np.cross(a,b)
        s = np.linalg.norm(v)
        c = np.dot(a,b.T)
        vx = np.matrix([[0,-v[0,2],v[0,1]],[v[0,2],0,-v[0,0]],[-v[0,1],v[0,0],0]])
        transform = np.eye(3,3) + vx + (vx*vx)*(1/(1+c[0,0]))
        return transform


    def center_and_align(self, parttype=1, align_percentile=10, return_transform=False, use_transform=False, trans=None, verbose=False, centeronly=False):
        """Center and align the particles in the region either with a supplied transformation matrix
            or by computation of the mean angular momentum of some range of particles (defined by some percentile)
            of the radii
        """
        ptypeind = self._ptypeind[parttype]
        if not self.centered:
            for i in range(len(self.coordinates)):
                self.coordinates[i] -= np.array(self.center)
        radii = np.linalg.norm(self.coordinates[ptypeind], axis=1)
        inside = radii < np.percentile(radii, align_percentile)
        self.bulkvel = np.median(self.velocities[ptypeind][inside], axis=0)
        for i in range(len(self.velocities)):
            self.velocities[i] -= np.array(self.bulkvel)
        if centeronly:
            return None
        self.centered = True
        if use_transform:
            t = trans
        else:
            t = self._transform(self.angular_momentum(parttype=parttype, percentile=align_percentile))
        if verbose:
            print('Transforming Coordinates...')
        for i in range(len(self.coordinates)):
            self.coordinates[i] = np.einsum('ij,aj->ai', t, self.coordinates[i])
            self.velocities[i] = np.einsum('ij,aj->ai', t, self.velocities[i])
        if return_transform:
            return t

    def _single_X_H(self,X,H,element):
        solar = self.solar_abundances[element]
        solarH = self.solar_abundances['Hydrogen']
        return np.log10(X/H)-np.log10(solar/solarH)


    def abundance_ratios(self,gas=False,smoothed=True):
        """ Compute element abundance ratios for the region, returns a dict of [X/H] """
        if smoothed:
            e_key = 'SmoothedElementAbundance'
        else:
            e_key = 'ElementAbundance'
        if gas:
            parttype = 0
        else:
            parttype = 4
        entries = []
        H = self.get_dataset(parttype,os.path.join(e_key,'Hydrogen'))
        for i in range(len(self.elements)):
            if self.elements[i] == 'Hydrogen' or self.elements[i] == 'Sulphur':
                continue
            X = self.get_dataset(parttype,os.path.join(e_key,self.elements[i]))
            entries.append((self.elements[i],self._single_X_H(X,H,self.elements[i])))
        return dict(entries)

class Subfind:
    """ Basic Subfind superclass which finds the relevant files.

    arguments:
        run - the run (e.g. L0012N0188)
        model - an EAGLE model (e.g. Ref)
        tag - a tag string specifying a snapshot output (e.g. 028_z000p000)

    history:
        written - Mackereth (UoB) - 22/11/2019
     """
    def __init__(self, run, model, tag):
        #store the snapshot identity info
        self.run = run
        self.model = model
        self.tag = tag
        if release == 'public':
            self.simlabel = self.model+self.run
            self.snaplabel = 'groups_'+self.tag
            self.base_subfile = 'eagle_subfind_tab_'+self.tag
            self.path = os.path.join(base_path, self.simlabel, self.snaplabel)
        elif release == 'ARI':
            self.snaplabel = 'groups_'+self.tag
            self.base_subfile = 'eagle_subfind_tab_'+self.tag
            self.path = os.path.join(base_path, self.run, self.model, 'data', self.snaplabel)
        else:
            raise Exception('private/custom data access is not yet implemented!')
        if not os.path.exists(os.path.join(self.path, self.base_subfile+'.0.hdf5')):
            raise Exception('could not see snapshot data in directory: '+self.path)
        #get the files related to this snapshot and load some of their metadata
        self.files = natural_sort(glob.glob(os.path.join(self.path, self.base_subfile+'*.hdf5')))
        self.nfiles = len(self.files)
        self.header_dict = dict(h5py.File(self.files[0], 'r')['/Header'].attrs.items())
        self.BoxSize = self.header_dict['BoxSize']
        self.HubbleParam  = self.header_dict['HubbleParam']
        self.Omega0, self.OmegaLambda, self.OmegaBaryon, self.a0 = self.header_dict['Omega0'], self.header_dict['OmegaLambda'], self.header_dict['OmegaBaryon'], self.header_dict['ExpansionFactor']
        self.datasets = {}
        bases = ['FOF', 'Subhalo']
        for base in bases:
            self.datasets[base] = list(h5py.File(self.files[0], 'r')[base].keys())

    def get_dataset(self, dataset, physical=False, cgs=False):
        """ get the data for a given entry in the HDF5 files """
        out = []
        if physical:
            #find conversion factor
            factor = self._conversion_factor(dataset, self.a0, self.HubbleParam, cgs=cgs)
        elif not physical and cgs:
            factor = h5py.File(self.files[0], 'r')[dataset].attrs['CGSConversionFactor']
        else:
            #else just multiply by 1!
            factor = 1
        for file in self.files:
            # load this file and get the particles
            out.append(np.array(h5py.File(file, 'r')[dataset])[:] * factor)
        return np.concatenate(out)


    def _conversion_factor(self, key, a, h, cgs=False):
        aexp_scale, h_scale = self._get_conversion_factor_exponents(key)
        if cgs:
            cgs_factor = h5py.File(self.files[0], 'r')[key].attrs['CGSConversionFactor']
        else:
            cgs_factor = 1.
        return a**(aexp_scale)*h**(h_scale)*cgs_factor

    def _get_conversion_factor_exponents(self, key):
        aexp_scale = h5py.File(self.files[0], 'r')[key].attrs['aexp-scale-exponent']
        h_scale = h5py.File(self.files[0], 'r')[key].attrs['h-scale-exponent']
        return aexp_scale, h_scale


def natural_sort(l):
    """natural sort using regex (adapted by Mark Byers on StackOverflow
    from http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html)"""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def _particle_type_present(type, file):
    head = dict(h5py.File(file, 'r')['/Header'].attrs.items())
    return head['NumPart_ThisFile'][type] > 0
