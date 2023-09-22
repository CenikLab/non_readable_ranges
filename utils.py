from ribopy.core.get_gadgets import get_region_boundaries, get_reference_names
from functools import cache
import numpy as np

@cache
def intevl(ribo_object, experiment_id):
    # ribo_object = ribo_data("%s.ribo"%experiment_id)
    data_tmp = ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data = data_tmp.iloc[6:26]
    pct_85 = sum(data["%s" % experiment_id]) * 0.85
    # pct_90=sum(data["%s"%experiment_id])*0.90
    value = data[data["%s" % experiment_id] == data["%s" % experiment_id].max()][
        "%s" % experiment_id
    ].values[0]
    mmin = mmax = data[data["%s" % experiment_id] == data["%s" % experiment_id].max()][
        "read_length"
    ].values[0]
    while value <= pct_85:
        if mmax < 40 and mmin > 21:
            if (
                data[data["read_length"] == mmax + 1]["%s" % experiment_id].values[0]
                >= data[data["read_length"] == mmin - 1]["%s" % experiment_id].values[0]
            ):
                mmax += 1
                value += data[data["read_length"] == mmax]["%s" % experiment_id].values[
                    0
                ]
            else:
                mmin -= 1
                value += data[data["read_length"] == mmin]["%s" % experiment_id].values[
                    0
                ]
        elif mmax == 40:
            mmin -= 1
            value += data[data["read_length"] == mmin]["%s" % experiment_id].values[0]
        elif mmin == 21:
            mmax += 1
            value += data[data["read_length"] == mmax]["%s" % experiment_id].values[0]
    # print(min,max)
    read_pct = value / sum(data["%s" % experiment_id])
    return mmin, mmax, read_pct

@cache
def get_cds_range_lookup(ribo):
    """
    Create a dict of gene to region ranges, so that the CDS range can be found for a given experiment.
    """
    names = get_reference_names(ribo._handle)
    if ribo.alias != None:
        names = map(ribo.alias.get_alias, names)
    boundaries = get_region_boundaries(ribo._handle)
    boundary_lookup = dict(zip(list(names), boundaries))
    return boundary_lookup


def cap_outliers(arr, thresh, filter_zeros=True):
    arr = arr.copy()
    if filter_zeros:
        arr = arr[arr > 0]
    if len(arr) == 0: return arr
    cap = np.percentile(arr, thresh)
    arr[arr > cap] = cap
    return arr

def cap_outliers_cds_only(arr, gene, boundary_lookup, thresh=99, filter_zeros=False):
    start, end = boundary_lookup[gene][1]
    arr = arr.copy()
    arr = arr[start:end]
    if filter_zeros:
        arr = arr[arr > 0]
    if len(arr) == 0: return arr
    cap = np.percentile(arr, thresh)
    arr[arr > cap] = cap
    return arr
