import eagleSqlTools as sql
import numpy as np
import os

username = os.environ['EAGLE_USERNAME']
password = os.environ['EAGLE_PASSWORD']

con = sql.connect(username, password=password)

model_dict = {'REFERENCE': 'Ref',
              'RECAL': 'Recal',
              'AGNdT9': 'AGNdT9'}

class SnapShotInfo:
    """ simple class to contain information from the database about snapshots

    history:
        written - Mackereth (UoB) - 05/02/2020
    """
    def __init__(self):
        query = 'SELECT\
        snap.Snapnum as snapnum,\
        snap.Redshift as redshift,\
        snap.LookbackTime as tlookback,\
        snap.LumDistance as lumdistance\
        FROM\
        Snapshots as snap'
        snapinfo = sql.execute_query(con,query)
        self.Snapnum = snapinfo['snapnum']
        self.Redshift = snapinfo['redshift']
        self.LookbackTime = snapinfo['tlookback']
        self.LumDistance = snapinfo['lumdistance']

class MergerTree:
    """ class to construct the merger tree of a given z=0 galaxy.

    arguments:
        ID - galaxy identifier, can be integer GalaxyID or tuple of [GroupNumber,SubGroupNumber,SnapNum]

    history:
        written - Mackereth (UoB) - 05/02/2020
    """
    def  __init__(self,ID, model='REFERENCE', run='L0100N1504'):
        self.model = model_dict[model]
        self.run = run
        if isinstance(ID, int):
            self.MainDescendantID = ID
            #query to get other descendant info:
            query = f'SELECT\
            gal.Snapnum as snapnum,\
            gal.GroupNumber as groupnum,\
            gal.SubGroupNumber as subgroupnum,\
            gal.TopLeafID as TopLeafID\
            FROM\
            {model_dict[model]+run}_Subhalo as gal\
            WHERE\
            gal.GalaxyID = {ID}'
            result = sql.execute_query(con,query)
            self.MainDescendantGroupNumber = result['groupnum']
            self.MainDescendantSubGroupNumber = result['subgroupnum']
            self.MainDescendantSnapNum = result['snapnum']
            self.MainDescendantTopLeafID = result['TopLeafID']

        else:
            # is [GroupNumber,SubGroupNumber,SnapNum]
            self.MainDescendantGroupNumber = ID[0]
            self.MainDescendantSubGroupNumber = ID[1]
            self.MainDescendantSnapNum = ID[2]
            query = f'SELECT\
            gal.GalaxyID as galaxyid,\
            gal.TopLeafID as TopLeafID\
            FROM\
            {model_dict[model]+run}_Subhalo as gal\
            WHERE\
            gal.GroupNumber = {ID[0]} and\
            gal.SubGroupNumber = {ID[1]} and\
            gal.SnapNum = {ID[2]}'
            result = sql.execute_query(con,query)
            self.MainDescendantID = result['galaxyid']
            self.MainDescendantTopLeafID = result['TopLeafID']

    def get_all_progenitors(self, mass_limit=1e9):
        """ return mass, coordinates and redshift of all progenitors """
        query = f'SELECT\
                 prog.GalaxyID as GalaxyID,\
                 prog.MassType_Star as mstar,\
                 prog.CentreOfPotential_x as x,\
                 prog.CentreOfPotential_y as y,\
                 prog.CentreOfPotential_z as z,\
                 prog.SnapNum as snapnum,\
                 prog.Redshift as redshift,\
                 prog.DescendantID as DescendantID,\
                 prog.LastProgID as LastProgID,\
                 prog.TopLeafID as  TopLeafID\
                 FROM\
                 {self.model+self.run}_Subhalo as des,\
                 {self.model+self.run}_Subhalo as prog\
                 WHERE\
                 des.GalaxyID = {self.MainDescendantID} and\
                 des.SnapNum = {self.MainDescendantSnapNum} and\
                 prog.MassType_star > {mass_limit} and\
                 prog.GalaxyID between {self.MainDescendantID} and des.LastProgID'
        result = sql.execute_query(con, query)
        #identify main branch
        mainbranch = result['TopLeafID'] == self.MainDescendantTopLeafID
        #store all progeny off the main branch (we want these separate for a few reasons...)
        self.ProgenitorIDs = result['GalaxyID'][~mainbranch]
        self.ProgenitorStellarMass = result['mstar'][~mainbranch]
        self.ProgenitorCoPs = np.dstack([result['x'], result['y'], result['z']])[0][~mainbranch]
        self.ProgenitorSnapNums = result['snapnum'][~mainbranch]
        self.ProgenitorRedshift = result['redshift'][~mainbranch]
        self.ProgenitorDescendantID = result['DescendantID'][~mainbranch]
        self.ProgenitorLastProgID = result['LastProgID'][~mainbranch]
        self.ProgenitorTopLeafID = result['TopLeafID'][~mainbranch]
        #store the Main Branch separately
        self.MainBranchIDs = result['GalaxyID'][mainbranch]
        self.MainBranchStellarMass = result['mstar'][mainbranch]
        self.MainBranchCoPs = np.dstack([result['x'], result['y'], result['z']])[0][mainbranch]
        self.MainBranchSnapNums = result['snapnum'][mainbranch]
        self.MainBranchRedshift = result['redshift'][mainbranch]
        self.MainBranchDescendantID = result['DescendantID'][mainbranch]
        self.MainBranchLastProgID = result['LastProgID'][mainbranch]
        self.MainBranchTopLeafID = result['TopLeafID'][mainbranch]
