import numpy as np

""" A great deal of this code is thieved from the Peano-Hilbert sorting algorithms included
in Gadget-2, written by Volker Springel"""

quadrants = [
  [[[0, 7], [1, 6]], [[3, 4], [2, 5]]],
  [[[7, 4], [6, 5]], [[0, 3], [1, 2]]],
  [[[4, 3], [5, 2]], [[7, 0], [6, 1]]],
  [[[3, 0], [2, 1]], [[4, 7], [5, 6]]],
  [[[1, 0], [6, 7]], [[2, 3], [5, 4]]],
  [[[0, 3], [7, 4]], [[1, 2], [6, 5]]],
  [[[3, 2], [4, 5]], [[0, 1], [7, 6]]],
  [[[2, 1], [5, 6]], [[3, 0], [4, 7]]],
  [[[6, 1], [7, 0]], [[5, 2], [4, 3]]],
  [[[1, 2], [0, 3]], [[6, 5], [7, 4]]],
  [[[2, 5], [3, 4]], [[1, 6], [0, 7]]],
  [[[5, 6], [4, 7]], [[2, 1], [3, 0]]],
  [[[7, 6], [0, 1]], [[4, 5], [3, 2]]],
  [[[6, 5], [1, 2]], [[7, 4], [0, 3]]],
  [[[5, 4], [2, 3]], [[6, 7], [1, 0]]],
  [[[4, 7], [3, 0]], [[5, 6], [2, 1]]],
  [[[6, 7], [5, 4]], [[1, 0], [2, 3]]],
  [[[7, 0], [4, 3]], [[6, 1], [5, 2]]],
  [[[0, 1], [3, 2]], [[7, 6], [4, 5]]],
  [[[1, 6], [2, 5]], [[0, 7], [3, 4]]],
  [[[2, 3], [1, 0]], [[5, 4], [6, 7]]],
  [[[3, 4], [0, 7]], [[2, 5], [1, 6]]],
  [[[4, 5], [7, 6]], [[3, 2], [0, 1]]],
  [[[5, 2], [6, 1]], [[4, 3], [7, 0]]]
]

rotxmap_table = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22]
rotymap_table = [1, 2, 3, 0, 16, 17, 18, 19,11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7]
rotx_table = [ 3, 0, 0, 2, 2, 0, 0, 1]
roty_table = [ 0, 1, 1, 2, 2, 3, 3, 0]
sense_table = [-1, -1, -1, 1, 1, -1, -1, -1]


def find_cell_coords(cell_size, x, y, z):
    #find the coordinates of the cell corresponding to position x,y,z
    nx, ny, nz = x//cell_size, y//cell_size, z//cell_size
    return nx.astype(int), ny.astype(int), nz.astype(int)

def ph_key(x,y,z,bits=6):
    #return the P-H key for the position x, y, z
    mask = 1 << (bits - 1)
    key = 0
    rotation = 0
    sense = 1
    i = 0
    while i < bits:
        if x & mask: bitx = 1
        else: bitx = 0
        if y & mask: bity = 1
        else: bity = 0
        if z & mask: bitz = 1
        else: bitz = 0
        quad = quadrants[rotation][bitx][bity][bitz]
        key = key << 3
        if sense == 1:
            key += quad
        else:
            key += 7-quad
        rotx = rotx_table[quad]
        roty = roty_table[quad]
        sense *= sense_table[quad]
        while rotx > 0:
            rotation = rotxmap_table[rotation]
            rotx -= 1
        while roty > 0:
            rotation = rotymap_table[rotation]
            roty -= 1
        i += 1
        mask = mask >> 1
    return key

def coordinate_grid(centre, side_length, n=70):
    #make a grid of coordinates in a cube around a given point
    xstart = centre[0]-1.*side_length/2.
    ystart = centre[1]-1.*side_length/2.
    zstart = centre[2]-1.*side_length/2.
    xend = centre[0]+1.*side_length/2.
    yend = centre[1]+1.*side_length/2.
    zend = centre[2]+1.*side_length/2.
    return np.mgrid[xstart:xend:n*1j, ystart:yend:n*1j, zstart:zend:n*1j].reshape(3,n*n*n).T

def get_unique_grid_keys(grid, cell_size, box_size, bits=6):
    #get the unique P-H keys corresponding to a grid of points
    grid = np.mod(grid, box_size) #wrap regions outside simulation volume
    nxyz = np.array(find_cell_coords(cell_size,grid[:,0], grid[:,1],grid[:,2]))
    nxyz = np.unique(nxyz.T, axis=0)
    keys = np.zeros(len(nxyz))
    for i in range(len(nxyz)):
        keys[i] = ph_key(nxyz[i,0], nxyz[i,1], nxyz[i,2], bits=bits)
    keys = np.unique(keys)
    return keys
