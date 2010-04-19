import numpy as np
import pyfits
import pysynphot as S
from pysynphot import observationmode
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

def getDate():
    """ Returns a formatted string with the current date.
        [This is simply a copy of 'getDate()' from pytools.fileutil]
    """
    _ltime = _time.localtime(_time.time())
    date_str = _time.strftime('%Y-%m-%dT%H:%M:%S',_ltime)

    return date_str

def createPrimaryHDU(filename,numpars,parnames,instrument):
    """ Create a Primary Header for the multi-extension FITS reference table
    """
    d = observationmode.getref()
    phdu = pyfits.PrimaryHDU()
    phdu.header.update('date',getDate(),comment="Date FITS file was generated")
    phdu.header.update('filename',filename,comment='name of file')
    phdu.header.update('nextend',3,comment='number of extensions in file')
    phdu.header.update('photzpt',-21.1,comment='Photometric zero-point for STMAG system')
    phdu.header.update('parnum',numpars,comment='Number of parameterized variables')
    for p in range(1,numpars+1):
        phdu.header.update('par'+str(p)+'name',parnames[p],comment='Name of 1st parameterized variable,if any')
    phdu.header.update('dbtable','IMPHTTAB')
    phdu.header.update('instrume',instrument)
    phdu.header.update('synswver',S.__version__,comment='Version of synthetic photometry software')
    phdu.header.update('graphtab',d['graphtable'],comment='HST Graph table')
    phdu.header.update('comptab',d['comptable'],comment='HST Components table')
    phdu.header.update('useafter','')
    phdu.header.update('pedigree','Test data')
    phdu.header.updata('descrip','photometry keywords reference file')

    return phdu
    
    
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
    obsmode_vals = []
    nmode_vals = []
    ped_vals =['Test data']*nfilt
    descrip_vals = ['Generated Apr 2010 for testing only']*nfilt
    
    filtname_vals = []
    fpars_vals = []
    npar_vals = []
    flam_datacol_vals = []
    plam_datacol_vals = []
    bw_datacol_vals = []
    for filt in filters:
        # For each filter combination (row in the table)...
        filtname,fpars = parseFilters(filt)
        # keep track of how many parameterized variables are used in this obsmode
        npars = len(fpars)
        npar_vals.append(npars)
        filtname_vals.append(filtname)
        fpars_vals.append(fpars)
        
        if npars == 0: nstr = ''
        else: nstr = str(npars)
        flam_datacol_vals.append('PHOTFLAM'+nstr)
        plam_datacol_vals.append('PHOTPLAM'+nstr)
        bw_datacol_vals.append('PHOTBW'+nstr)

        
    # Start by determining the maximum number of parameters in any given obsmode
    max_npars = np.array(npar_vals,np.int32).max()
    nelem_vals = [[]]*max_npars+1
    parvals_vals = [[]]*max_npars+1
    for filtname,fpars,npars in zip(filtname_vals,fpars_vals,npar_vals):
        #extract entries from 'filtdata' for only the values given in 'fpars'
        for f,n in zip(fpars,range(1,max_npars+1)):
            nelem = len(filtdata[f])
            nelem_vals[n].append(nelem)
            parvals_vals[n].append(filtdata[f])
    
        # insure that every row has the same number of entries for each column
        nelem_vals[0].append(0)
        parvals_vals[0].append(None)
        for i in range(n,max_npars):
            nelem_vals[i].append(0)
            parvals_vals[i].append([0.0])
        
    flam_vals = [[]]*max_npars+1
    plam_vals = [[]]*max_npars+1
    bw_vals = [[]]*max_npars+1
    for filtname,fpars,npars in zip(filtname_vals,fpars_vals,npar_vals):
        # define obsmode specific dictionary of parameterized variables
        fdict = {}
        for f in fpars:
            fdict[f] = filtdata[f]

        # Now build up list of all obsmodes with all combinations of 
        # parameterized variables values
        omode = basemode+','+filtname
        olist = generateObsmodes(omode,fdict)
        # interpret first obsmode for use as OBSMODE column value in table
        obsmode_vals.append(interpretObsmode(olist[0]))        
        
        # Use these obsmodes to generate all the values needed for the row
        nmodes = len(olist)
        nmode_vals.append(nmodes)
        photflam = np.zeros(nmodes,np.float64)
        photplam = np.zeros(nmodes,np.float64)
        photbw = np.zeros(nmodes,np.float64)
        
        for obsmode,n in zip(olist,range(nmodes)):
            value = computeValues(obsmode)
            photflam[n] = value['PHOTFLAM']
            photplam[n] = value['PHOTPLAM']
            photbw[n] = value['PHOTBW']
        for o in zip(range(max_npars+1)):
            if o == npars:
                flam_vals[o].append(photflam)
                plam_vals[o].append(photplam)
                bw_vals[o].append(photbw)
            else:
                flam_vals[o].append(0.0)
                plam_vals[o].append(0.0)
                bw_vals[o].append(0.0)
        del photflam,photplam,photbw
        
    # Finally, create the structures needed to define this row in the FITS table
    nelem_arr = np.array(nelem_vals,np.int32)
    # compute the maximum number of values in any results column
    nmodes_max = np.array(nmode_vals,np.int32).max()
    
    # Define each column in the table based on max_npars which are not different
    # from one extension to the other
    obsmode_col = Column(name='obsmode',format='40A',array=np.array(obsmode_vals))
    pedigree_col = Column(name='pedigree',format='30A',array=np.array(ped_vals))
    descrip_col = Column(name='descrip',format='67A',array=np.array(descrip_vals))
    datacol_col = {}
    datacol_col['PHOTFLAM'] = Column(name='datacol',format='12A',array=np.array(flam_datacol_vals))
    datacol_col['PHOTPLAM'] = Column(name='datacol',format='12A',array=np.array(plam_datacol_vals))
    datacol_col['PHOTBW'] = Column(name='datacol',format='12A',array=np.array(bw_datacol_vals))
    
    parvals_cols = []
    nelem_cols = []
    # for each parameterized element, create a set of columns specifying the
    # range of values for that parameter and the number of elements covering that range
    # namely, the PAR<n>VALUES and NELEM<n> columns
    for p in range(1,max_npars+1):
        nelem_cols.append(Column(name="NELEM"+str(p),format="I",array=np.array(nelem_vals[p],np.int16)))
        parvals_cols.append(Column(name="PAR"+str(p)+"VALUES",format="PD[]",array=np.array(parvals_vals[p],np.float64)))
        
    # create the set of results columns
    flam_cols = []
    plam_cols = []
    bw_cols = []
    for p in range(max_npars+1):
        if p == 0:
            format_str = 'D'
            pstr = ''
        else:
            format_str = 'PD[]'
            pstr = str(p)
        flam_cols.append(Column(name='PHOTFLAM'+pstr,format=format_str,array=np.array(flam_vals[p],np.float64)))
        plam_cols.append(Column(name='PHOTPLAM'+pstr,format=format_str,array=np.array(plam_vals[p],np.float64)))
        bw_cols.append(Column(name='PHOTBW'+pstr,format=format_str,array=np.array(bw_vals[p],np.float64)))

    # Now create the FITS file with the table in each extension
    phdu = createPrimaryHDU(output,2,['MJD','FR'],'acs')
            
    ftab = pyfits.HDUList()
    ftab.append(phdu)
    ftab.append(pyfits.new_table([obsmode_col,datacol_col['PHOTFLAM']]+flam_cols+parvals_cols+nelem_cols+[pedigree_col,descrip_col]))
    ftab.append(pyfits.new_table([obsmode_col,datacol_col['PHOTPLAM']]+plam_cols+parvals_cols+nelem_cols+[pedigree_col,descrip_col]))
    ftab.append(pyfits.new_table([obsmode_col,datacol_col['PHOTBW']]+bw_cols+parvals_cols+nelem_cols+[pedigree_col,descrip_col]))
    ftab.writeto(output)
    