import numpy as np

def kappa_co(coordinates, velocities, mass, rlim_kpc=30):
    """
    kappa_co following Correa et al. (2017)

    arguments:
        coordinates - coordinates in Mpc
        velocities - velocities in km s^{-1}
        mass - stellar mass (1e10 M_sun)

    history:
        written - Mackereth (UoB) - 22/07/2020
    """
    r = np.linalg.norm(coordinates, axis=1)*1e3
    mask = r < 30
    mass = snap.get_dataset(4,'Mass')[mask]
    pos = coordinates[mask]*1e3
    vel = velocities[mask]

    Rxy = np.linalg.norm(pos[:,:2], axis=1)
    L = (mass*np.cross(pos, vel).T).T
    K = np.sum(0.5*mass*np.linalg.norm(vel,axis=1)**2)
    Krot = np.sum(0.5*mass*(L[:,2]/(mass*Rxy))**2)
    return Krot/K

def t_con(dens, m_g=1.8e6):
    mu = 1.22
    m_p = 1.6726e-24
    T_0 = 8e3
    xH = 0.752
    rho_0 = 0.1 * m_p / xH
    P_0_on_kb = rho_0*T_0/mu/m_p
    P_eos_on_kb = P_0_on_kb*((dens*6.77e-31)/rho_0)**(4./3.)
    P_eos = P_eos_on_kb*1.38e-16
    G = 6.67e-8 #cgs
    A = 1.515e-4 / 1e6
    nhs = dens2nh(dens)
    ns = np.zeros(len(dens))
    t_g = 1.67e9*(P_eos_on_kb/10**3)**-0.2 #schaye and dalla vecchia 08 approximation
    return t_g, P_eos

def dens2nh(dens):
    tdens = (0.752/1.6726e-24)*(dens*6.77e-31)
    return tdens


def angular_momentum(coordinates, velocities, mass):
    """
    Compute the angular momentum of some coordinates and velocities
    """
    vec = np.cross(coordinates,velocities*mass[:,np.newaxis])
    tot = np.sum(vec, axis=0)
    return tot/np.linalg.norm(tot)

def _transform(vector):
    """Build a transformation matrix"""
    a = vector
    b = np.matrix([0,0,1])
    v = np.cross(a,b)
    s = np.linalg.norm(v)
    c = np.dot(a,b.T)
    vx = np.matrix([[0,-v[0,2],v[0,1]],[v[0,2],0,-v[0,0]],[-v[0,1],v[0,0],0]])
    transform = np.eye(3,3) + vx + (vx*vx)*(1/(1+c[0,0]))
    return transform


def center_and_align(coordinates,velocities,mass, center=None, transform=None, align_percentile=10, verbose=False, return_transform=False):
    """Center and align the particles in the region either with a supplied transformation matrix
        or by computation of the mean angular momentum of some range of particles (defined by some percentile)
        of the radii
    """
    if center is None:
        center = find_center_of_mass(coordinates)
    tcoordinates = coordinates - np.array(center)
    radii = np.linalg.norm(tcoordinates, axis=1)
    inside = radii < np.percentile(radii, align_percentile)
    bulkvel = np.median(velocities[inside], axis=0)
    tvelocities = velocities  - np.array(bulkvel)
    if transform is None:
        transform = _transform(angular_momentum(tcoordinates, tvelocities,mass))
    if verbose:
        print('Transforming Coordinates...')
    tcoordinates = np.einsum('ij,aj->ai', transform, tcoordinates)
    tvelocities = np.einsum('ij,aj->ai', transform, tvelocities)
    if return_transform:
        return tcoordinates, tvelocities, transform
    else:
        return tcoordinates, tvelocities

def find_center_of_mass(coords, quantile=5, max_iter=20, final_npart=1000):
    """ find the center of mass of a set of particles """
    percentiles = np.percentile(coords, [quantile,100-quantile], axis=0)
    mask = np.all(coords > percentiles[0], axis=1) & np.all(coords < percentiles[1], axis=1)
    i = 1
    while sum(mask) > final_npart and i < max_iter:
        coords = coords[mask]
        percentiles = np.percentile(coords, [quantile,100-quantile], axis=0)
        mask = np.all(coords > percentiles[0], axis=1) & np.all(coords < percentiles[1], axis=1)
        i += 1
    return np.median(coords, axis=0)
