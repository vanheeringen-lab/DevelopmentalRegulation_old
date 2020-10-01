"""
All utility functions for the analysis which we do not want to see in the
jupyter notebooks
"""
import itertools

import numpy as np


# google spreadsheet url of samples
samples_url = "https://docs.google.com/spreadsheets/d/1X2YiJ-ZlX5usBGFr2L_-Jmt_NkLDJeu8m2jq4OjyEXI/export?gid=664155010&format=tsv"

# ass2name and name2ass are dictionaries that convert assembly 
# names to human readable names and vice versa
def ass2name(assembly):
    """
    assembly name to human readable
    """
    return _ass2name[assembly]

def name2ass(name):
    """
    human readable name to assembly name
    """
    return {v: k for k, v in _ass2name.items()}[name]

_ass2name = {
    "Astyanax_mexicanus-2.0": "cavefish",
    "Bter_1.0": "Bumblebee",
    "ARS-UCD1.2": "cow",
    "BraLan2": "amphioxus",
    "Bl71nemr": "amphioxus",
    "ce11": "nematode",
    "KH": "sea vase",
    "danRer11": "zebrafish",
    "dm6": "fruit fly",
    "GRCg6a": "chicken",
    "hydra": "hydra",
    "GRCm39": "mouse",
    "mm10": "mouse",
    "ASM20922v1": "sea anemone",
    "ASM223467v1": "medaka",
    "Phmamm_MTP2014": "sea squirt",
    "HLtupMer6": "tegu",
    "ASM318616v1": "turbot",
    "Spur_5.0": "sea urchin",
    "xenTro9": "frog",
}


def _local_maximum(in_array, radius=3):
    """
    private function to determine the maxima in a 2d array
    """
    diam = 2 * radius + 1
    pad_array = np.pad(in_array, radius, 'minimum')
    
    # threshold = max(0, np.quantile(in_array, 0.0))
    threshold = 0
    out_array = in_array > threshold
    for x_off, y_off in itertools.product(range(diam), repeat=2):
        x_stop = x_off-diam + 1 if x_off-diam + 1 != 0 else None
        y_stop = y_off-diam + 1 if y_off-diam + 1 != 0 else None
        out_array &= in_array >= pad_array[slice(x_off, x_stop), slice(y_off, y_stop)]
    return out_array


def get_outline(array, radius=4):
    """
    Get the outline for ... 
    """
    maxima = _local_maximum(array, radius)
    
    ver_seg = np.where(maxima[:,1:] != maxima[:,:-1])
    hor_seg = np.where(maxima[1:,:] != maxima[:-1,:])

    l = []
    for p in zip(*hor_seg):
        l.append((p[1], p[0]+1))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan,np.nan))

    # and the same for vertical segments
    for p in zip(*ver_seg):
        l.append((p[1]+1, p[0]))
        l.append((p[1]+1, p[0]+1))
        l.append((np.nan, np.nan))

    # now we transform the list into a numpy array of Nx2 shape
    segments = np.array(l)

    # now we need to know something about the image which is shown
    #   at this point let's assume it has extents (x0, y0)..(x1,y1) on the axis
    #   drawn with origin='lower'
    # with this information we can rescale our points
    offset = -0.5
    x0 = 0 + offset
    #x1 = len(df1.columns) + offset
    x1 = array.shape[0] + offset
    y0 = 0 + offset
    y1 = array.shape[1] + offset

    segments[:,0] = x0 + (x1-x0) * segments[:,0] / (x1 - x0)
    segments[:,1] = y0 + (y1-y0) * segments[:,1] / (y1 - y0)
    
    return segments[:,0], segments[:,1], maxima
