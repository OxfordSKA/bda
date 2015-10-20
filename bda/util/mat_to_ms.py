#!/usr/bin/env python

# Run using:
# casapy --nologger --nogui --log2term -c mat_to_ms <MAT name> <input MS template> <output MS name>

import os
import sys
import scipy
import scipy.io
import collections
import numpy as np

# http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries
def loadmat(filename):
    '''
    this function should be called instead of direct scipy.io.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=False)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], np.ndarray):
            for i in dict[key]:
                if isinstance(i[0], scipy.io.matlab.mio5_params.mat_struct):
                    dict[key] = _todict(i[0])
    return dict

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = collections.OrderedDict()
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        #print strg, type(elem), elem.shape
        if isinstance(elem, np.ndarray) and elem.size > 0:
            if isinstance(elem[0][0], scipy.io.matlab.mio5_params.mat_struct):
                dict[strg] = _todict(elem[0][0])
            else:
                dict[strg] = elem
        elif isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def add_model_data():
    # Add the model data column and initialise to zero.
    model_desc = {
        'MODEL_DATA': {
            'comment': 'model data',
            'dataManagerGroup': 'ModelTiled',
            'dataManagerType': 'TiledShapeStMan',
            'maxlen': 0,
            'ndim': 2,
            'option': 0,
            'valueType': 'complex'
        }
    }
    model_dminfo = {
        '*6': {
            'COLUMNS': np.array(['MODEL_DATA'], dtype='|S11'),
            'NAME': 'ModelTiled',
            'SEQNR': 5,
            'SPEC': {
                'ActualMaxCacheSize': 0,
                'DEFAULTTILESHAPE': np.array([1, 1, 64262], dtype=np.int32),
                'HYPERCUBES': {
                    '*1': {
                        'BucketSize': 514096,
                        'CellShape': np.array([1, 1], dtype=np.int32),
                        'CubeShape': np.array([1, 1, 1606550], dtype=np.int32),
                        'ID': {},
                        'TileShape': np.array([1, 1, 64262], dtype=np.int32)
                    }
                },
                'IndexSize': 1,
                'MAXIMUMCACHESIZE': 0,
                'SEQNR': 5
            },
            'TYPE': 'TiledShapeStMan'
        }
    }
    tb.addcols(model_desc, model_dminfo)
    tb.putcol('MODEL_DATA', np.zeros((1, 1, tb.nrows()), dtype='c16'))

def add_corrected_data():
    # Add the corrected data column and initialise to zero.
    cor_desc = {
        'CORRECTED_DATA': {
            'comment': 'corrected data',
            'dataManagerGroup': 'CorrectedTiled',
            'dataManagerType': 'TiledShapeStMan',
            'maxlen': 0,
            'ndim': 2,
            'option': 0,
            'valueType': 'complex'
            }
        }
    cor_dminfo = {
        '*7': {
            'COLUMNS': np.array(['CORRECTED_DATA'], dtype='|S15'),
            'NAME': 'CorrectedTiled',
            'SEQNR': 6,
            'SPEC': {
                'ActualMaxCacheSize': 0,
                'DEFAULTTILESHAPE': np.array([1, 1, 64262], dtype=np.int32),
                'HYPERCUBES': {
                    '*1': {
                        'BucketSize': 514096,
                        'CellShape': np.array([1, 1], dtype=np.int32),
                        'CubeShape': np.array([1, 1, 1606550], dtype=np.int32),
                        'ID': {},
                        'TileShape': np.array([1, 1, 64262], dtype=np.int32)
                    }
                },
                'IndexSize': 1,
                'MAXIMUMCACHESIZE': 0,
                'SEQNR': 6
            },
            'TYPE': 'TiledShapeStMan'
        }
    }
    tb.addcols(cor_desc, cor_dminfo)
    tb.putcol('CORRECTED_DATA', np.zeros((1, 1, tb.nrows()), dtype='c16'))

def copy_ms(ms_in, ms_out):
    """Make a copy of a MS without copying any rows of the main table."""
    tb.open(ms_in, nomodify=True)
    tbcopy = tb.copy(ms_out, deep=True, valuecopy=True, norows=True,
                     returnobject=False)
    if tbcopy:
        tbcopy.close()

    # Copy all subtables intact.
    sub_tables = tb.getkeywords()
    tb.close()
    for s in sorted(sub_tables):
        if str(sub_tables[s]).startswith('Table'):
            if s == 'SORTED_TABLE':
                continue
            print 'Copying', ms_in + '/' + s, 'to', ms_out + '/' + s
            tb.open(ms_in + '/' + s)
            tbcopy = tb.copy(ms_out + '/' + s, deep=True, valuecopy=True,
                             norows=False, returnobject=False)
            if tbcopy:
                tbcopy.close()
            tb.close()

if __name__ == '__main__':
    # Get MAT file name and MS names from command line.
    mat_name = sys.argv[-3]
    template_name = sys.argv[-2]
    ms_name = sys.argv[-1]

    # Load the MAT file into a dictionary.
    d = loadmat(mat_name)

    # Copy the input (template) MS to a new one.
    copy_ms(template_name, ms_name)

    # Open the output MS.
    print "Opening table", ms_name
    tb.open(ms_name, nomodify=False)

    # Check if corrected data or model data columns should be added.
    if 'CORRECTED_DATA' in d and not 'CORRECTED_DATA' in tb.colnames():
        print "Adding corrected data"
        add_corrected_data()
    if 'MODEL_DATA' in d and not 'MODEL_DATA' in tb.colnames():
        print "Adding model data"
        add_model_data()

    # Iterate the dictionary (all the columns).
    found_length = 0
    for k in d:
        # Check if column exists in the main table.
        if k in tb.colnames():
            req_length = d[k].shape[-1]
            if found_length == 0:
                if req_length == 0:
                    continue
                print "Adding rows", req_length
                tb.addrows(req_length)
                found_length = 1
            elif req_length > 0 and req_length != tb.nrows():
                raise RuntimeError("Inconsistent row dimension!")
            if req_length > 0:
                print "Putting column", k, type(d[k]), d[k].shape, d[k].dtype
                tb.putcol(k, d[k])

    print "Closing table"
    tb.close()
