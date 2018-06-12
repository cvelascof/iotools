"""Methods for reading files containing 2d precipitation fields.

The methods in this module implement the following interface:

  read_xxx(filename, optional arguments)

where xxx is the name (or abbreviation) of the file format and filename is the 
name of the input file.

The output of each method is a three-element tuple containing the two-dimensional 
precipitation field and georeferencing and metadata dictionaries.

The geodata dictionary contains the following mandatory key-value pairs:
  ll_lon       longitude of the lower-left corner of the data raster
  ll_lat       latitude of the lower-left corner of the data raster
  ur_lon       longitude of the upper-right corner of the data raster
  ur_lat       latitude of the upper-right corner of the data raster
  projection   PROJ.4-compatible projection definition
  xpixelsize   grid resolution in x-direction (meters)
  ypixelsize   grid resolution in y-direction (meters)
  yorigin      a string specifying the origin of the y-axis:
               'upper' = upper border
               'lower' = lower border

The metadata dictionary contains the following mandatory key-value pairs:
  missingval   value to indicate missing data
"""

import gzip
from matplotlib.pyplot import imread
from numpy import nan
import pyproj

def read_pgm(filename, dtype=float, gzipped=False, convert_dbz_to_r=True, 
             dbz_to_r_a=223.0, dbz_to_r_b=1.53):
    """Read a 8-bit PGM radar reflectivity composite from the FMI archive and 
    optionally convert the reflectivity values to precipitation rates.

    Parameters
    ----------
    filename : str
        Name of the file to read from.
    dtype : type
        The output datatype for the dataset that is read from the file.
    gzipped : bool
        If True, the input file is treated as a compressed gzip file.
    convert_dbz_to_r : bool
        If True, the reflectivities Z are converted to precipitation rates R via 
        the formula Z=a*R^b.
    dbz_to_r_a : float
        Multiplier for the Z->R conversion.
    dbz_to_r_b : float
        Exponent for the Z->R conversion.
  
    Returns
    -------
    out : tuple
        A two-element tuple containing the precipitation field read from the PGM 
        file and the associated georeferencing data.
    """
    metadata = _read_pgm_metadata(filename, gzipped=gzipped)

    if gzipped == False:
        R = imread(filename)
    else:
        R = imread(gzip.open(filename, 'r'))
    geodata = _read_pgm_geodata(metadata)

    MASK = R == metadata["missingval"]
    R = R.astype(dtype)
    R[MASK] = nan
    R = (R - 64.0) / 2.0
    if convert_dbz_to_r:
        R = pow(pow(10.0, R / 10.0) / dbz_to_r_a, 1.0 / dbz_to_r_b)

    return R, geodata, metadata

def _read_pgm_geodata(metadata):
    geodata = {}

    projdef = ""

    if metadata["type"][0] != "stereographic":
        raise ValueError("unknown projection %s" % metadata["type"][0])
    projdef += "+proj=stere "
    projdef += " +lon_0=" + metadata["centrallongitude"][0] + 'E'
    projdef += " +lat_0=" + metadata["centrallatitude"][0] + 'N'
    projdef += " +lat_ts=" + metadata["truelatitude"][0]
    # These are hard-coded because the projection definition is missing from the 
    # PGM files.
    projdef += " +a=6371288"
    projdef += " +x_0=380886.310"
    projdef += " +y_0=3395677.920"
    projdef += " +no_defs"
    #
    geodata["projection"] = projdef
  
    ll_lon,ll_lat = [float(v) for v in metadata["bottomleft"]]
    ur_lon,ur_lat = [float(v) for v in metadata["topright"]]

    pr = pyproj.Proj(projdef)
    x1,y1 = pr(ll_lon, ll_lat)
    x2,y2 = pr(ur_lon, ur_lat)

    geodata["ll_lon"] = ll_lon
    geodata["ll_lat"] = ll_lat
    geodata["ur_lon"] = ur_lon
    geodata["ur_lat"] = ur_lat

    geodata["x1"] = x1
    geodata["y1"] = y1
    geodata["x2"] = x2
    geodata["y2"] = y2

    geodata["xpixelsize"] = float(metadata["metersperpixel_x"][0])
    geodata["ypixelsize"] = float(metadata["metersperpixel_y"][0])

    geodata["yorigin"] = "upper"
  
    return geodata

def _read_pgm_metadata(filename, gzipped=False):
    metadata = {}
  
    if gzipped == False:
        f = open(filename, 'r')
    else:
        f = gzip.open(filename, 'r')
  
    l = f.readline()
    while l[0] != '#':
        l = f.readline()
    while l[0] == '#':
        x = l[1:].strip().split(' ')
        if len(x) >= 2:
            k = x[0]
            v = x[1:]
            metadata[k] = v
        else:
            l = f.readline()
            continue
        l = f.readline()
    l = f.readline()
    metadata["missingval"] = int(l)
    f.close()
  
    return metadata
