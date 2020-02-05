import eagleSqlTools as sql

username = os.environ['EAGLE_USERNAME']
password = os.environ['EAGLE_PASSWORD']

con = sql.connect(username, password=password)

model_dict = {'REFERENCE': 'Ref',
              'RECAL': 'Recal',
              'AGNdT9': 'AGNdT9'}


class SnapShotInfo:
    """ simple class to contain information from the database about snapshots

    history:
        written - Mackereth (UoB) - 22/11/2019
     """
     def __init__(self):
         query = '''SELECT
            snap.Snapnum as snapnum,
            snap.Redshift as redshift,
            snap.LookbackTime as tlookback,
            snap.LumDistance as lumdistance
           FROM
            Snapshots as snap'''
        snapinfo = sql.execute_query(con,query)
        self.Snapnum = snapinfo['snapnum']
        self.Redshift = snapinfo['redshift']
        self.LookbackTime = snapinfo['tlookback']
        self.LumDistance = snapinfo['lumdistance']
