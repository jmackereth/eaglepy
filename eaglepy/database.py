import eagleSqlTools as sql
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
        if isinstance(ID, int):
            self.DescendantID = ID
            #query to get other descendant info:
            query = 'SELECT\
            gal.Snapnum as snapnum,\
            gal.GroupNumber as groupnum,\
            gal.SubGroupNumber as subgroupnum\
            FROM\
            %s%s_Subhalo as gal\
            WHERE\
            gal.GalaxyID = %s' % (model_dict[model], run, ID)
            result = sql.execute_query(con,query)
            self.DescendantGroupNumber = result['groupnum'][0]
            self.DescendantSubGroupNumber = result['subgroupnum'][0]
            self.DescendantSnapNum = result['snapnum'][0]
        else:
            # is [GroupNumber,SubGroupNumber,SnapNum]
            self.DescendantGroupNumber = ID[0]
            self.DescendantSubGroupNumber = ID[1]
            self.DescendantSnapNum = ID[2]
            query = 'SELECT\
            gal.GalaxyID as galaxyid\
            FROM\
            %s%s_Subhalo as gal\
            WHERE\
            gal.GroupNumber = %s and\
            gal.SubGroupNumber = %s and\
            gal.SnapNum = %s' % (model_dict[model], run, ID[0], ID[1], ID[2])
            result = sql.execute_query(con,query)
            self.DescendantID = result['galaxyid'][0]
