#!/usr/bin/python

import os
import shutil
from numpy import array, int32
import time

def create_main_table(cal_table, ms):
    desc = {
        'TIME': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'double'
        },
        'FIELD_ID': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'int'
        },
        'SPECTRAL_WINDOW_ID': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'int'
        },
        'ANTENNA1': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'int'
        },
        'ANTENNA2': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'int'
        },
        'INTERVAL': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'double'
        },
        'SCAN_NUMBER': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'int'
        },
        'OBSERVATION_ID': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'option': 5,
            'valueType': 'int'
        },
        'CPARAM': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'ndim': -1,
            'option': 0,
            'valueType': 'complex'
        },
        'PARAMERR': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'ndim': -1,
            'option': 0,
            'valueType': 'float'
        },
        'FLAG': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'ndim': -1,
            'option': 0,
            'valueType': 'boolean'
        },
        'SNR': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'ndim': -1,
            'option': 0,
            'valueType': 'float'
        },
        'WEIGHT': {
            'comment': '',
            'dataManagerGroup': 'MSMTAB',
            'dataManagerType': 'StandardStMan',
            'maxlen': 0,
            'ndim': -1,
            'option': 0,
            'valueType': 'float'
        },
        '_define_hypercolumn_': {}
    }
    dminfo = {
        '*1':
        {
            'COLUMNS': array(['ANTENNA1', 'ANTENNA2', 'CPARAM', 'FIELD_ID',
                              'FLAG', 'INTERVAL', 'OBSERVATION_ID',
                              'PARAMERR', 'SCAN_NUMBER', 'SNR',
                              'SPECTRAL_WINDOW_ID', 'TIME', 'WEIGHT'],
                              dtype='|S19'),
            'NAME': 'MSMTAB',
            'SEQNR': 0,
            'SPEC': {'ActualCacheSize': 2,
                     'BUCKETSIZE': 2560,
                     'IndexLength': 9974,
                     'PERSCACHESIZE': 2},
            'TYPE': 'StandardStMan'
        }
    }
    tb.create(tablename=cal_table, tabledesc=desc, dminfo=dminfo, nrow=197*200)

    path = os.path.realpath(cal_table)
    keywords = {
        'ParType': 'Complex',
        'MSName': '%s' % ms,
        'VisCal': 'G Jones',
        'PolBasis': 'unknown',
        'OBSERVATION': 'Table: %s/OBSERVATION' % path,
        'ANTENNA': 'Table: %s/ANTENNA' % path,
        'FIELD': 'Table: %s/FIELD' % path,
        'SPECTRAL_WINDOW': 'Table: %s/SPECTRAL_WINDOW' % path,
        'HISTORY': 'Table: %s/HISTORY' % path
    }
    # for key in keywords:
    #     print key, keywords[key]
    #     tb.putkeyword(key, keywords[key])
    tb.putkeyword(keyword='ANTENNA',
                  value='Table: TMP.cal/ANTENNA')

    tb.flush()
    tb.close()

def create_antenna_table(cal_table):
    desc = {'DISH_DIAMETER': {'comment': 'Physical diameter of dish',
                       'dataManagerGroup': 'MSMTAB',
                       'dataManagerType': 'StandardStMan',
                       'maxlen': 0,
                       'option': 0,
                       'valueType': 'double'},
     'FLAG_ROW': {'comment': 'Flag for this row',
                  'dataManagerGroup': 'MSMTAB',
                  'dataManagerType': 'StandardStMan',
                  'maxlen': 0,
                  'option': 0,
                  'valueType': 'boolean'},
     'MOUNT': {'comment': 'Mount type e.g. alt-az, equatorial, etc.',
               'dataManagerGroup': 'MSMTAB',
               'dataManagerType': 'StandardStMan',
               'maxlen': 0,
               'option': 0,
               'valueType': 'string'},
     'NAME': {'comment': 'Antenna name, e.g. VLA22, CA03',
              'dataManagerGroup': 'MSMTAB',
              'dataManagerType': 'StandardStMan',
              'maxlen': 0,
              'option': 0,
              'valueType': 'string'},
     'OFFSET': {'comment': 'Axes offset of mount to FEED REFERENCE point',
                'dataManagerGroup': 'MSMTAB',
                'dataManagerType': 'StandardStMan',
                'maxlen': 0,
                'ndim': 1,
                'option': 5,
                'shape': array([3], dtype=int32),
                'valueType': 'double'},
     'POSITION': {'comment': 'Antenna X,Y,Z phase reference position',
                  'dataManagerGroup': 'MSMTAB',
                  'dataManagerType': 'StandardStMan',
                  'maxlen': 0,
                  'ndim': 1,
                  'option': 5,
                  'shape': array([3], dtype=int32),
                  'valueType': 'double'},
     'STATION': {'comment': 'Station (antenna pad) name',
                 'dataManagerGroup': 'MSMTAB',
                 'dataManagerType': 'StandardStMan',
                 'maxlen': 0,
                 'option': 0,
                 'valueType': 'string'},
     'TYPE': {'comment': 'Antenna type (e.g. SPACE-BASED)',
              'dataManagerGroup': 'MSMTAB',
              'dataManagerType': 'StandardStMan',
              'maxlen': 0,
              'option': 0,
              'valueType': 'string'},
     '_define_hypercolumn_': {}}

    dminfo = {'*1': {'COLUMNS': array(['DISH_DIAMETER', 'FLAG_ROW', 'MOUNT', 'NAME', 'OFFSET', 'POSITION',
       'STATION', 'TYPE'],
      dtype='|S14'),
        'NAME': 'MSMTAB',
        'SEQNR': 0,
        'SPEC': {'ActualCacheSize': 2,
                 'BUCKETSIZE': 3332,
                 'IndexLength': 174,
                 'PERSCACHESIZE': 2},
        'TYPE': 'StandardStMan'}}


    tb.create(tablename=cal_table+'/ANTENNA', tabledesc=desc, dminfo=dminfo,
              nrow=197)
    tb.close()



def main():
    # --------------------------------------------------------
    cal_table = 'TMP.cal'
    ms = 'test_cor.ms'
    # --------------------------------------------------------


    if os.path.isdir(cal_table):
        print 'Removing old table'
        shutil.rmtree(cal_table)

    t0 = time.time()
    create_main_table(cal_table, ms)
    create_antenna_table(cal_table)
    print time.time()-t0


if __name__ == "__main__":
    main()
