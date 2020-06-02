#function to delete the visual properties produced by EDC2 on each run in order to decluster the viewport
# import functions used within this function
from abaqus import *
from abaqusConstants import *

import part

def axisplanedel(model, part) :
    # import model and part information
    p = mdb.models[model].parts[part]
    # initialise counters
    k = len(p.features)
    i = 1
    l = 1
    #Delete all planes from the viewport
    while l < k:
        if 'plane' in p.features.keys()[i]:
            key1 = p.features.keys()[i]
            del p.features[key1]
        else:
            i += 1
        l += 1
    # import model and part information
    p = mdb.models[model].parts[part]
    # initialise counters
    k = len(p.features)
    i = 1
    l = 1
    #Delete all axis from the viewport
    while l < k:
        if 'axis' in p.features.keys()[i]:
            key1 = p.features.keys()[i]
            del p.features[key1]
        else:
            i += 1
        l += 1