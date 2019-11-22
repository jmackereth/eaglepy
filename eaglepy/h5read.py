import h5py
import os
import glob
import re
import numpy as np
from . import peano

base_path = os.environ['EAGLE_BASE_PATH']
release = os.environ['EAGLE_ACCESS_TYPE']

class Snapshot:
    def __init__(self, run, model, tag):
        #store the snapshot identity info
        self.run = run
        self.model = model
        self.tag = tag
        if release == 'public':
            self.simlabel = self.model+self.run
            self.snaplabel = 'snapshot_'+self.tag
            self.base_subfile = 'snap_'+self.tag
            self.path = os.path.join(base_path, self.simlabel, self.snaplabel)
        else:
            raise Exception('private/custom data access is not yet implemented!')
        if not os.path.exists(os.path.join(self.path, self.base_subfile+'.0.hdf5')):
            raise Exception('could not see snapshot data in directory: '+self.path)
        #get the files related to this snapshot and load some of their metadata
        self.files = natural_sort(glob.glob(os.path.join(self.path, self.base_subfile+'*')))
        self.nfiles = len(self.files)
        self.header_dict = dict(h5py.File(self.files[0])['/Header'].attrs.items())
        self.BoxSize = self.header_dict['BoxSize']
        self.HubbleParam  = self.header_dict['HubbleParam']
        self.NumPartTotal = self.header_dict['NumPart_Total']
        self.ParticleTypePresent = np.where(self.NumPartTotal > 0)[0]
        self._ptypeind = {self.ParticleTypePresent[i]:i for i in range(len(self.ParticleTypePresent))}
        #get the Hash Table info for P-H key sorting
        self.HashBits = dict(h5py.File(self.files[0])['/HashTable'].attrs.items())['HashBits']
        self.HashGridSideLength = 2**self.HashBits
        self.HashGridCellSize = self.BoxSize/self.HashGridSideLength
        self.firstkeys = np.zeros((len(self.ParticleTypePresent),self.nfiles))
        self.lastkeys  = np.zeros((len(self.ParticleTypePresent),self.nfiles))
        self.datasets = {}
        for ii,parttype in enumerate(self.ParticleTypePresent):
            self.firstkeys[ii] = np.array(h5py.File(self.files[0])['/HashTable/PartType'+str(parttype)+'/FirstKeyInFile'])
            self.lastkeys[ii] = np.array(h5py.File(self.files[0])['/HashTable/PartType'+str(parttype)+'/LastKeyInFile'])
            self.datasets['PartType'+str(parttype)] = list(h5py.File(self.files[0])['/PartType'+str(parttype)].keys())


class SnapshotRegion(Snapshot):
    def __init__(self, run, model, tag, center, sidelength):
        super().__init__(run, model, tag)
        self.center = center
        self.sidelength = sidelength
        self._index_region(self.center, self.sidelength)

    def _index_region(self, center, side_length, phgrid_n=70):
        """ Load a region defined by a central cordinate and a side length
        arguments:
        center - the [x,y,z] coordinate of the desired center (simulation units)
        side_length - the desired side length (in the simulation units)

        keyword arguments:
        phgrid_n - the number of grid points along a side length to look for PH cells (default 70)
        """
        #work out which files contain the desired region
        grid = peano.coordinate_grid(center, side_length, n=phgrid_n)
        keys  = peano.get_unique_grid_keys(grid, self.HashGridCellSize, self.BoxSize, bits=self.HashBits)
        self.files_for_region = []
        coordinates = []
        indices = []
        for ii,type in enumerate(self.ParticleTypePresent):
            Nfiles = self._get_parttype_files(type, keys)
            self.files_for_region.append(np.array(self.files)[Nfiles])
            #now load the coordinates in these files and save the indices for each particle type
            thistypecoord, thistypeindices = self._get_parttype_indices(type, self.files_for_region[ii])
            coordinates.append(thistypecoord)
            indices.append(thistypeindices)
        self.coordinates = coordinates
        self.indices = indices


    def _get_parttype_indices(self, parttype, files):
        """get the coordinates and indices for a given particle type in a given region"""
        coords, indices = [], []
        for file in files:
            # load the file
            thisfilecoords = np.array(h5py.File(file)['/PartType'+str(parttype)+'/Coordinates'])
            # mask it to the region desired
            mask = (np.fabs(thisfilecoords[:,0]-self.center[0]) < self.sidelength/2.) &\
                    (np.fabs(thisfilecoords[:,1]-self.center[1]) < self.sidelength/2.) &\
                     (np.fabs(thisfilecoords[:,2]-self.center[2]) < self.sidelength/2.)
            #store the coordinates and the indices of these particles in the file
            thisfileindices = np.where(mask)[0]
            coords.append(thisfilecoords[mask])
            indices.append(thisfileindices)
        return np.concatenate(coords), indices

    def get_dataset(self, parttype, dataset):
        """ get the data for a given entry in the HDF5 file for the given region """
        key = os.path.join('/PartType'+str(parttype),dataset)
        out = []
        ptypeind = self._ptypeind[parttype]
        for ii,file in enumerate(self.files_for_region[ptypeind]):
            # load this file and get the particles
            out.append(np.array(h5py.File(file)[key])[self.indices[ptypeind][ii]])
        return np.concatenate(out)

    def _get_parttype_files(self, parttype, keys):
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





def natural_sort(l):
    #natural sort using regex (adapted by Mark Byers on StackOverflow
    #from http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html)
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)
