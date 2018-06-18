"""Utilities for finding archived files that match the given criteria."""

from datetime import datetime, timedelta
import fnmatch
import os


def find_by_date(date, root_path, path_fmt, fn_pattern, fn_ext, timestep, 
                 num_prev_files=0):
    """List input files whose timestamp matches the given date.

    Parameters
    ----------
    date : datetime.datetime
        The given date.
    root_path : str
        The root path to search the input files.
    path_fmt : str
        Path format. It may consist of directory names separated by '/' and 
        date/time specifiers beginning with '%' (e.g. %Y/%m/%d).
    fn_pattern : str
        The name pattern of the input files without extension. The pattern can 
        contain time specifiers (e.g. %H, %M and %S).
    fn_ext : str
        Extension of the input files.
    timestep : int
        Time step between consecutive input files (minutes).
    num_prev_files : int
        Optional, number of previous files to find before the given timestamp.
  
    Returns
    -------
    out : tuple
        If num_prev_files=0, return a pair containing the found file name and the 
        corresponding timestamp as a datetime.datetime object. Otherwise, return 
        a tuple of two lists, the first one for the file names and the second one 
        for the correspondign timestemps. The lists are sorted in ascending order 
        with respect to timestamp.
    """
    filenames  = []
    timestamps = []
  
    for i in xrange(num_prev_files+1):
        curdate = date - timedelta(minutes=i*timestep)
        fn = _find_matching_filename(curdate, root_path, path_fmt, fn_pattern, fn_ext)
        filenames.append(fn)

        timestamps.append(curdate)
  
    if num_prev_files > 0:
        return (filenames[::-1], timestamps[::-1])
    else:
        return (filenames[0], timestamps[0])

def _find_matching_filename(date, root_path, path_fmt, fn_pattern, fn_ext):
    path = _generate_path(date, root_path, path_fmt)
    fn = None
  
    if os.path.exists(path):
        fn = datetime.strftime(date, fn_pattern) + '.' + fn_ext
        
        # test for wildcars
        if '?' in fn:
            filenames = os.listdir(path)
            if len(filenames) > 0:
                for filename in filenames:
                    if fnmatch.fnmatch(filename, fn):
                        fn = filename
                        break
            
        fn = os.path.join(path, fn)
        fn = fn if os.path.exists(fn) else None
  
    return fn

def _generate_path(date, root_path, path_fmt):
    f = lambda t: datetime.strftime(date, t) if t[0] == '%' else t
    if path_fmt != "":
        tokens = [f(t) for t in path_fmt.split('/')]
        subpath = os.path.join(*tokens)
    
        return os.path.join(root_path, subpath)
    else:
        return root_path
