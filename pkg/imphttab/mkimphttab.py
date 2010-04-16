import numpy as np
import pyfits
import pysynphot as S
from pyfits import Column

def computeValues(obsmode):
    """ Compute the 3 photometric values needed for a given obsmode string
        using pysynphot
        This routine will return a dictionary with the photometry keywords and
        the computed value of the keyword
    """    
    # Define the bandpass for this obsmode
    bp = S.ObsBandpass(obsmode)
    # set up the flat spectrum used for the computing the photometry keywords
    sp = S.FlatSpectrum(1,fluxunits='flam')
    # create the observation using these elements
    obs = S.Observation(sp,bp)
    # compute the photometric values
    valdict = {}
    
    valdict['PHOTFLAM'] = obs.effstim('flam')/obs.effstim('counts')
    valdict['PHOTPLAM'] = obs.pivot()
    valdict['PHOTBW'] = bp.rmswidth()
    
    return valdict
    
def generateObsmodes(basemode, pardict):
    """ Generate a set of obsmode strings spanning all the combinations of the 
        parameter as specified in the input dictionary
        
        The input will be a dictionary with the name of each parameterized variabble
        as the key and a list of values for that key; eg., 
            {'mjd':[52334,53919.99,53920,55516],'fr853n':[8158.0,8531.5,8905.0]}

        The basemode will be the root of the obsmode string including all 
        non-parameterized components; eg., acs,wfc1,f850lp
    """
    obsmode_str = '%s,%s#%0.4f'
    olist = []
    if len(pardict) == 0:
        # we don't have any parameterized variables, so just return the basemode
        olist.append(basemode)
    elif len(pardict) == 1:
        # build up list of obsmodes covering all combinations of paramterized
        # variable values
        key = pardict.keys()[0]
        for val in pardict[key]:
            olist.append(obsmode_str%(basemode,key,val))
    else:        
        nkeys = len(pardict)
        for nkey in range(nkeys-1):
            key = pardict.keys()[nkey]
            for val in pardict[key]:
                pdict = {}
                for k in pardict.keys()[nkey+1:]:
                    pdict[k] = pardict[k]
                olist.extend(generateObsmodes(obsmode_str%(basemode,key,val),pdict))

    return olist

def interpretObsmode(obsmode):
    ''' Convert a full obsmode string with parameterized values into a string
        which only lists the parameterized variable names without values

        The return value corresponds to the OBSMODE value for each row in the
        table.
        
        For example, convert:
            acs,wfc1,mjd#52334.0000,fr853n#8158.0000
        into:
            acs,wfc1,mjd,fr853n
    '''
    ospl = obsmode.split(',')
    omode = ''
    for o in ospl:
        if '#' in o: o = o.split('#')[0]
        omode += o+','
    omode = omode.rstrip(',')
    return omode

        
def read_dict(fname):
    '''read a python dictionary from a file that was written
       with write_dict. [COPIED directly from pyetc/etc_engine/util.py]
    '''
    f=open(fname,'r')
    datastr = f.read()
    f.close()
    # convert DOS file to Unix - otherwise the eval will fail
    datastr = datastr.replace('\r','')
    try :
        datadict = eval(datastr)
    except Exception, e:
        print 'EXCEPTION:',e
        print 'cannot eval data in file ',fname
        raise
    return datadict
    
def parseFilters(filters):
    """ Parse the filters specification and return a string with any 
        non-parameterized filter names, and a list of parameterized filter names.
    """
    fspl = filters.split(',')
    fpars = []
    fval = ''
    for f in fspl:
        if f.find('#') > -1:
            f = f.replace('#','')
            fpars.append(f)
        fval += '%s,'%f
    fval.rstrip(',')
    
    return fval,fpars

def buildIMPHTTAB(output,basemode,filters,filtdata,clobber=True):
    """ Create a (sample) IMPHTTAB file for a specified base 
        configuration (basemode) and a set of filter combinations (filters).
        
        The input 'filters' will be a list of strings, with parameterized 
        filters ending in '#', for all the combinations of filters to be 
        computed for the table with one row per combination.  
        
        The range of values spanned by each of the parameterized filters 
        will be specified in the external file or dictionary 'filtdata'.
        If 'filtdata' is the name of a file, the data must be formatted as
        a Python list of dictionaries, one dictionary per filter.
    """
    # check status of output file
    if os.path.exists(output):
        if clobber: os.remove(output)
        else: raise IOError,'Output file already exists. Please delete/rename before restarting.'
    # interpret input data
    # The 'filtdata' dict contains the values for ALL the parameterized variables
    # used in the obsmodes supported by this table
    if isinstance(filtdata,str):
        filtdata = read_dict(filtdata)
    
    # start building obsmodes for each row
    nfilt = len(filters)
    flam_vals = []
    plam_vals = []
    bw_vals = []
    npar_vals = []
    obsmode_vals = []
    ped_vals =['Test data']*nfilt
    descrip_vals = ['Generated Apr 2010 for testing only']*nfilt
    nelem_vals = []
    parvals_vals = []
    
    for filt in filters:
        # For each filter combination (row in the table)...
        filtname,fpars = parseFilters(filt)
        # keep track of how many parameterized variables are used in this obsmode
        npars = len(fpars)
        npar_vals.append(npars)
        
        #extract entries from 'filtdata' for only the values given in 'fpars'
        fdict = {}
        nelem = []
        parvals = []
        for f in fpars:
            fdict[f] = filtdata[f]
            nelem.append(len(fdict[f]))
            parvals.append(fdict[f])
        nelem_vals.append(nelem)
        parvals_vals.append(parvals)
        
        # Now build up list of all obsmodes with all combinations of 
        # parameterized variables values
        omode = basemode+','+filtname
        olist = generateObsmodes(omode,fdict)
        # interpret first obsmode for use as OBSMODE column value in table
        obsmode_vals.append(interpretObsmode(olist[0]))        
        
        # Use these obsmodes to generate all the values needed for the row
        nmodes = len(olist)
        photflam = np.zeros(nmodes,np.float64)
        photplam = np.zeros(nmodes,np.float64)
        photbw = np.zeros(nmodes,np.float64)
        
        for obsmode,n in zip(olist,range(nmodes)):
            value = computeValues(obsmode)
            photflam[n] = value['PHOTFLAM']
            photplam[n] = value['PHOTPLAM']
            photbw[n] = value['PHOTBW']
        flam_vals.append(photflam)
        plam_vals.append(photplam)
        bw_vals.append(photbw)
        del photflam,photplam,photbw
        
    # Finally, create the structures needed to define this row in the FITS table

    # Start by determining the maximum number of parameters in any given obsmode
    max_npars = np.array(npar_col,np.int32).max()

    nelem_arr = np.array(nelem_vals,np.int32)
    
    # Define each column in the table based on max_npars which are not different
    # from one extension to the other
    obsmode_col = Column(name='obsmode',format='40A',array=np.array(obsmode_vals))
    pedigree_col = Column(name='pedigree',format='30A',array=np.array(ped_vals))
    descrip_col = Column(name='descrip',format='67A',array=np.array(descrip_vals))

    parvals_cols = []
    for p in range(nfilt):
        pass

        