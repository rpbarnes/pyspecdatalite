import paramset
from os import listdir,environ
if paramset.myparams['figlist_type'] == 'figlistl':
    environ['ETS_TOOLKIT'] = 'qt4'
    import matplotlib; matplotlib.use('Agg')
import textwrap
import matplotlib.transforms as mtransforms
from numpy import sqrt as np_sqrt
from numpy.lib.recfunctions import rename_fields,drop_fields
from mpl_toolkits.mplot3d import axes3d
from matplotlib.collections import PolyCollection
from matplotlib.colors import LightSource
from matplotlib.lines import Line2D
from scipy.interpolate import griddata as scipy_griddata
import tables
import warnings
from inspect import ismethod
from numpy.core import rec
from matplotlib.pyplot import cm
from copy import deepcopy 
import traceback
import sympy
from scipy.optimize import leastsq
from scipy.signal import fftconvolve
import scipy.sparse as sparse
import numpy.lib.recfunctions as recf
from inspect import getargspec
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from pylab import *
#rc('image',aspect='auto',interpolation='bilinear') # don't use this, because it gives weird figures in the pdf
rc('image',aspect='auto',interpolation='nearest')
rcParams['xtick.direction'] = 'out'
rcParams['xtick.major.size'] = 12
rcParams['xtick.minor.size'] = 6
rcParams['ytick.direction'] = 'out'
rcParams['ytick.major.size'] = 12
rcParams['ytick.minor.size'] = 6
#rcParams['lines.linewidth'] = 3.0
rcParams['legend.fontsize'] = 12
rcParams['axes.grid'] = False
rcParams['font.size'] = 18
#{{{ constants
k_B = 1.380648813e-23
mu_0 = 4e-7*pi
mu_B = 9.27400968e-24#Bohr magneton
epsilon_0 = 8.854187817e-12
hbar = 6.6260695729e-34/2./pi
N_A = 6.02214179e23
gamma_H = 4.258e7
#}}}

def mybasicfunction(first_figure = None):
    r'''this gives the format for doing the image thing
note also
nextfigure(fl,'name')
and
nextfigure({'lplotproperty':value})
'''
    fl = figlistini_old(first_figure)
    return figlistret(first_figure,figurelist,other)

def issympy(x):
    'tests if something is sympy (based on the module name)'
    return repr(x.__class__.__module__)[1:-1].split('.')[0] == 'sympy'

#{{{ function trickery
def mydiff(data,axis = -1):
    '''this will replace diff with a version that has the same number of indices, with the last being the copy of the first'''
    newdata = zeros(shape(data),dtype = data.dtype)
    indices = [slice(None,None,None)]*len(data.shape)
    indices[axis] = slice(None,-1,None)
    newdata[indices] = diff(data,axis = axis)
    #setfrom = list(indices)
    #indices[axis] = -1
    #setfrom[axis] = 0
    #newdata[indices] = newdata[setfrom]
    return newdata
#}}}
def normal_attrs(obj):
    myattrs = filter(lambda x: not ismethod(obj.__getattribute__(x)),dir(obj))
    myattrs = filter(lambda x: not x[0:2] == '__',myattrs)
    return myattrs
def showtype(x):
    if type(x) is ndarray:
        return ndarray,x.dtype
    else:
        return type(x)
def emptyfunction():
    pass

def lookup_rec(A,B,indexpair):
    r'''look up information about A in table B (i.e. chemical by index, etc)
    indexpair is either the name of the index
    or -- if it's differently named -- the pair of indices
    given in (A,B) respectively
    
    This will just drop any fields in B that are also in A,
    and the output uses the first indexname
    
    note that it it seems like the join_rec function above may be more efficient!!'''
    raise RuntimeError('You should now use decorate_rec!!')
    if type(indexpair) not in [tuple,list]:
        indexpair = (indexpair,indexpair)
    Bini = copy(B)
    B = recf.drop_fields(B,( set(B.dtype.names) & set(A.dtype.names) ) - set([indexpair[1]])) # indexpair for B gets dropped later anyways
    joined = []
    for j in A:
        matchedrows =  B[B[indexpair[1]] == j[indexpair[0]]]
        for matchedrow in matchedrows:
            joined.append((j,matchedrow))
    if len(joined) == 0:
        raise CustomError('Unable to find any matches between',A[indexpair[0]],'and',B[indexpair[1]],'!')
    whichisindex = joined[0][1].dtype.names.index(indexpair[1])
    allbutindex = lambda x: list(x)[0:whichisindex]+list(x)[whichisindex+1:]
    joined = concatenate([array(tuple(list(j[0])+allbutindex(j[1])),
                    dtype = dtype(j[0].dtype.descr+allbutindex(j[1].dtype.descr))).reshape(1) for j in joined])
    return joined
def reorder_rec(myarray,listofnames,first = True):
    try:
        indices_to_move = [myarray.dtype.names.index(j) for j in listofnames]
    except:
        stuff_not_found = [j for j in listofnames if j not in myarray.dtype.names]
        if len(stuff_not_found) > 0:
            raise CustomError(stuff_not_found,'is/are in the list you passed, but not one of the fields, which are',myarray.dtype.names)
        else:
            raise CustomError('unknown problem')
    old_type = list(myarray.dtype.descr)
    new_type = [old_type[j] for j in indices_to_move] + [old_type[j] for j in range(0,len(old_type)) if j not in indices_to_move]
    new_list_of_data = [myarray[j[0]] for j in new_type]
    return rec.fromarrays(new_list_of_data,dtype = new_type)

def lambda_rec(myarray,myname,myfunction,*varargs):
    r'''make a new field "myname" which consists of "myfunction" evaluated with the fields given by "myargs" as arguments
    the new field is always placed after the last argument name
    if myname is in myargs, the original row is popped'''
    if len(varargs) == 1:
        myargs = varargs[0]
    elif len(varargs) == 0:
        myargs = [myname]
    else:
        raise CustomError("For the fourth argument, you must pass either a list with the names of the arguments, or nothing (to use the field itself as an argument)")
    if type(myargs) is str:
        myargs = (myargs,)
    if type(myargs) is not tuple:
        myargs = tuple(myargs)
    argdata = map((lambda x: myarray[x]),myargs)
    try:
        newrow = myfunction(*tuple(argdata))
    except TypeError:
        newrow = array([myfunction(*tuple([x[rownumber] for x in argdata])) for rownumber in range(0,len(argdata[0]))])
    if type(newrow) is list and type(newrow[0]) is str:
        newrow = array(newrow,dtype = '|S100')
    try:
        new_field_type = list(newrow.dtype.descr[0])
    except AttributeError:
        raise CustomError("evaluated function on",argdata,"and got back",newrow,"which appears not to be a numpy array")
    new_field_type[0] = myname
    starting_names = myarray.dtype.names
    #{{{ make the dtype
    new_dtype = list(myarray.dtype.descr)
    #{{{ determine if I need to pop one of the existing rows due to a name conflict
    eliminate = None
    if myname in myargs:
        eliminate = myname
        insert_before = starting_names.index(myname) # if we are replacing, we want it in the same place
        new_dtype = [j for j in new_dtype if j[0] != eliminate]
    #}}}
    # if we haven't already eliminated, determine where to put it
    if eliminate is None:
        insert_before = starting_names.index(myargs[-1])+1
    # insert the new field where it goes
    new_dtype.insert(insert_before,tuple(new_field_type))
    #}}}
    #{{{ separate starting_names and ending_names
    if eliminate is None:
        ending_names = starting_names[insert_before:]
        starting_names = starting_names[:insert_before]
    else: # if I'm eliminating, I don't want to include the eliminated one
        ending_names = starting_names[insert_before+1:]
        starting_names = starting_names[:insert_before]
    #}}}
    return rec.fromarrays([myarray[x] for x in starting_names if x != eliminate]+[newrow]+[myarray[x] for x in ending_names if x != eliminate],dtype = new_dtype)

def join_rec((A,a_ind),(B,b_ind)):
    raise RuntimeError('You should now use decorate_rec!!')

def decorate_rec((A,a_ind),(B,b_ind),drop_rows = False,verbose = False):
    r'''Decorate the rows in A with information in B --> if names overlap,
    keep the ones in A
    b_ind and a_ind can be either a single key, or a list of keys;
    if more than one element in B matches that in A, include both options!!'''
    #A = A.copy() # because I play with it later
    dropped_rows = None
    # first find the list of indices that give us the data we want
    #{{{ process the arguments
    if (type(b_ind) is str) and (type(a_ind) is str):
        b_ind = [b_ind]
        a_ind = [a_ind]
    if ((type(b_ind) is list) and (type(a_ind) is list)) and (len(b_ind) == len(a_ind)):
        pass
    else:
        raise ValueError('If you call a list for b_ind and/or a_ind, they must match in length!!!')
    if any([x not in B.dtype.names for x in b_ind]):
        problem_index = [x for x in b_ind if x not in B.dtype.names]
        raise ValueError(repr(problem_index)+' not in second argument, which has fields'+repr(B.dtype.names)+'!!!')
    if any([x not in A.dtype.names for x in a_ind]):
        problem_index = [x for x in a_ind if x not in A.dtype.names]
        raise ValueError(repr(problem_index)+' not in first argument, which has fields'+repr(A.dtype.names)+'!!!')
    #}}}
    B_reduced = B[b_ind] # a version of B reduced to only include the keys
    B_reduced = reorder_rec(B_reduced,b_ind)# again, because it doesn't do this just based on the indexing
    A_reduced = A[a_ind] # same for A
    A_reduced = reorder_rec(A_reduced,a_ind)# again, because it doesn't do this just based on the indexing
    # now, I need to generate a mapping from the b_ind to a_ind
    field_mapping = dict(zip(b_ind,a_ind))
    # now I change the names so they match and I can compare them
    B_reduced.dtype.names = tuple([field_mapping[x] for x in B_reduced.dtype.names])
    #{{{ now find the list of indices for B that match each value of A
    old_B_reduced_names,old_B_reduced_types = tuple(zip(*tuple(B_reduced.dtype.descr)))
    B_reduced.dtype = dtype(zip(A_reduced.dtype.names,old_B_reduced_types))
    if A_reduced.dtype != B_reduced.dtype:
        B_reduced.dtype = dtype(zip(old_B_reduced_names,old_B_reduced_types))
        raise CustomError('The datatype of A_reduced=',A_reduced.dtype,'and B_reduced=',B_reduced.dtype,'are not the same, which is going to create problems!')
    try:
        list_of_matching = [nonzero(B_reduced == j)[0] for j in A_reduced]
    except:
        raise CustomError('When trying to decorate, A_reduced=',A_reduced,'with B_reduced=',B_reduced,'one or more of the following is an empty tuple, which is wrong!:',[nonzero(B_reduced == j) for j in A_reduced])
    if verbose: print "(decorate\\_rec):: original list of matching",list_of_matching
    length_of_matching = array([len(j) for j in list_of_matching])
    if verbose: print "(decorate\\_rec):: length of matching is",length_of_matching
    if any(length_of_matching == 0):
        if drop_rows:
            if drop_rows == 'return':
                dropped_rows = A[length_of_matching == 0].copy()
            else:
                dropped_rows = A_reduced[length_of_matching == 0]
                print r'{\color{red}Warning! decorate\_rec dropped fields in the first argument',lsafen(repr(zip(A_reduced.dtype.names * len(dropped_rows),dropped_rows.tolist()))),r'}'
            #{{{ now, remove all trace of the dropped fields
            A = A[length_of_matching != 0]
            list_of_matching = [j for j in list_of_matching if len(j)>0]
            length_of_matching = [len(j) for j in list_of_matching]
            #}}}
        else:
            raise CustomError('There is no data in the second argument that has',b_ind,'fields to match the',a_ind,'fields of the first argument for the following records:',A_reduced[length_of_matching == 0],'if this is correct, you can set the drop_rows = True keyword argument to drop these fields')
    # now, do a neat trick of stackoverflow to collapse a nested list
    # this gives just the indices in B that match the values of A
    list_of_matching = [j for i in list_of_matching for j in i]
    #}}}
    if verbose: print "(decorate\\_rec):: list of matching is",list_of_matching
    # now grab the data for these rows
    add_data = B[list_of_matching]
    #{{{ finally, smoosh the two sets of data together
    #{{{ Now, I need to replicate the rows that have multiple matchesjk
    if any(length_of_matching > 1):
        index_replication_vector = [k for j in range(0,len(length_of_matching))
                for k in [j]*length_of_matching[j]]
        retval = A[index_replication_vector]
    else:
        retval = A.copy()
    #}}}
    #{{{ add the new fields
    new_dtypes = [j for j in B.dtype.descr if j[0] not in A.dtype.names]
    if verbose: print "(decorate\\_rec):: new dtypes:",repr(new_dtypes)
    try:
        retval = newcol_rec(retval,new_dtypes)
    except:
        raise CustomError("Problem trying to add new columns with the dtypes",new_dtypes)
    #}}}
    if verbose: print "(decorate\\_rec):: add data:",repr(add_data)
    for name in dtype(new_dtypes).names:
        if verbose: print "(decorate\\_rec):: trying to add data for",name,':',add_data[name][:]
        retval[name][:] = add_data[name][:]
    #}}}
    if drop_rows == 'return':
        return retval,dropped_rows
    else:
        return retval

def newcol_rec(A,new_dtypes):
    r'''add new, empty (i.e. random numbers) fields to A, as given by new_dtypes
    --> note that there are deeply nested numpy functions to do this, but the options are confusing, and I think the way these work is efficient'''
    if type(new_dtypes) is dtype:
        new_dtypes = new_dtypes.descr
    elif type(new_dtypes) is tuple:
        new_dtypes = [new_dtypes]
    elif type(new_dtypes) is list:
        if type(new_dtypes[0]) is not tuple:
            new_dtypes = [tuple(new_dtypes)]
    retval = empty(A.shape,dtype = A.dtype.descr + new_dtypes)
    for name in A.dtype.names:
        retval[name][:] = A[name][:]
    return retval

def applyto_rec(myfunc,myarray,mylist,verbose = False):
    r'apply myfunc to myarray with the intention of collapsing it to a smaller number of values'
    if type(mylist) is not list and type(mylist) is str:
        mylist = [mylist]
    combined = []
    j = 0
    #{{{ make the list "combined", which I later concatenate
    while len(myarray) > 0:
        thisitem = myarray[0] # always grab the first row of what's left
        #{{{ initialize the empty new row
        if j == 0:
            newrow = thisitem.reshape(1)
        newrow = newrow.copy()
        #}}}
        #{{{ make a mask for all items that are identified as the same data
        # and copy the identical data to newrow
        mask = myarray[mylist[0]] == thisitem[mylist[0]]
        newrow[mylist[0]] = thisitem[mylist[0]]
        for k in range(1,len(mylist)):
            mask &= myarray[mylist[k]] == thisitem[mylist[k]]
            newrow[mylist[k]] = thisitem[mylist[k]]
        #}}}
        if verbose: print lsafen('(applyto rec): for row %d, I select these:'%j)
        myarray_subset = myarray[mask]
        if verbose: print lsafen('(applyto rec): ',repr(myarray_subset))
        other_fields = set(mylist)^set(thisitem.dtype.names)
        if verbose: print lsafen('(applyto rec): other fields are:',other_fields)
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = myfunc(myarray_subset[thisfield])
            except:
                raise CustomError("error in applyto_rec:  You usually get this when one of the fields that you have NOT passed in the second argument is a string.  The fields and types are:",repr(myarray_subset.dtype.descr))
        if verbose: print lsafen("(applyto rec): for row %d, I get this as a result:"%j,newrow)
        combined.append(newrow) # add this row to the list
        myarray = myarray[~mask] # mask out everything I have used from the original matrix
        if verbose: print lsafen("(applyto rec): the array is now",repr(myarray))
        j += 1
    #}}}
    combined = concatenate(combined)
    if verbose: print lsafen("(applyto rec): final result",repr(combined),"has length",len(combined))
    return combined

def meanstd_rec(myarray,mylist,verbose = False,standard_error = False):
    r'this is something like applyto_rec, except that it applies the mean and creates new rows for the "error," where it puts the standard deviation'
    if type(mylist) is not list and type(mylist) is str:
        mylist = [mylist]
    combined = []
    other_fields = set(mylist)^set(myarray.dtype.names)
    if verbose: print '(meanstd_rec): other fields are',lsafen(other_fields)
    newrow_dtype = [[j,('%s_ERROR'%j[0],)+j[1:]] if j[0] in other_fields else [j] for j in myarray.dtype.descr]
    newrow_dtype = [k for j in newrow_dtype for k in j]
    if verbose: print lsafen('(meanstd rec): other fields are:',other_fields)
    #{{{ make the list "combined", which I later concatenate
    j = 0
    while len(myarray) > 0:
        thisitem = myarray[0] # always grab the first row of what's left
        #{{{ initialize the empty new row
        newrow = zeros(1,dtype = newrow_dtype)
        #}}}
        #{{{ make a mask for all items that are identified as the same data
        # and copy the identical data to newrow
        mask = myarray[mylist[0]] == thisitem[mylist[0]]
        newrow[mylist[0]] = thisitem[mylist[0]]
        for k in range(1,len(mylist)):
            mask &= myarray[mylist[k]] == thisitem[mylist[k]]
            newrow[mylist[k]] = thisitem[mylist[k]]
        #}}}
        if verbose: print lsafen('(meanstd rec): for row %d, I select these:'%j)
        myarray_subset = myarray[mask]
        if verbose: print lsafen('(meanstd rec): ',repr(myarray_subset))
        for thisfield in list(other_fields):
            try:
                newrow[thisfield] = mean(myarray_subset[thisfield])
                if standard_error:
                    newrow[thisfield+"_ERROR"] = std(myarray_subset[thisfield])/sqrt(len(myarray_subset[thisfield]))
                else:
                    newrow[thisfield+"_ERROR"] = std(myarray_subset[thisfield])
            except:
                raise CustomError("error in meanstd_rec:  You usually get this when one of the fields that you have NOT passed in the second argument is a string.  The fields and types are:",repr(myarray_subset.dtype.descr))
            #print 'for field',lsafe(thisfield),'I find',lsafen(newrow[thisfield])
        if verbose: print lsafen("(meanstd rec): for row %d, I get this as a result:"%j,newrow)
        combined.append(newrow) # add this row to the list
        myarray = myarray[~mask] # mask out everything I have used from the original matrix
        if verbose: print lsafen("(meanstd rec): the array is now",repr(myarray))
        j += 1
    #}}}
    combined = concatenate(combined)
    if verbose: print lsafen("(meanstd rec): final result",repr(combined),"has length",len(combined))
    return combined

def make_rec(*args,**kwargs):
    r'input,names or a single argument, which is a dictionary\nstrlen = 100 gives length of the strings (which need to be specified in record arrays)\nyou can also specify (especially useful with the dictionary format) the list order = [str1,str2,...] which orders the output records with the field containing str1 first, then the field containing str2, then any remaining fields'
    strlen = 100
    if 'strlen' in kwargs.keys():
        strlen = kwargs.pop('strlen')
    if 'order' in kwargs.keys():
        order = kwargs.pop('order')
    else:
        order = None
    if 'zeros_like' in kwargs.keys():
        zeros_like = kwargs.pop('zeros_like')
    else:
        zeros_like = False
    if len(kwargs)>0:
        raise CustomError("You have kwargs I don't understand!:",kwargs)
    if len(args) == 1 and (type(args[0]) is dict):
        names = args[0].keys()
        input = args[0].values()
    elif len(args) == 2:
        input = args[0]
        names = args[1]
    else:
        raise CustomError("I don't understand the arguments you passed to make_rec!!!\nshould be (list of values, list of field names), or a dictionary")
    #{{{ apply the order kwarg
    if order is not None:
        newindices = []
        for orderitem in order:
            newindices += [j for j,k in enumerate(names) if (k.find(orderitem)>-1 and j not in newindices)]
        newindices += [j for j,k in enumerate(names) if j not in newindices]
        names = [names[j] for j in newindices]
        input = [input[j] for j in newindices]
    #}}}
    if not (type(input) is list and type(names) is list):
        raise CustomError('you must enter a list for both')
    types = map(type,input)
    shapes = map(shape,input)
    if all([j == shapes[0] for j in shapes]):
        if shapes[0] == ():# if it's one dimensional
            equal_shapes = False
            shapes = [(1)]*len(shapes)
        else:
            equal_shapes = True
            shape_of_array = shapes[0]
            shapes = [()]*len(shapes)
    else:
        equal_shapes = False
    for j,k in enumerate(input):
        if type(k) is list and equal_shapes:
            k = k[0]
        if type(k) is str:
            types[j] = '|S%d'%strlen
        if type(k) is ndarray:
            types[j] = k.dtype
    try:
        mydtype = dtype(zip(names,types,shapes))
    except:
        raise CustomError('problem trying to make names',names,' types',types,'shapes',shapes)
    if zeros_like:
        retval = zeros(zeros_like,dtype = mydtype)
        return retval
    if equal_shapes:
        retval = empty(shape_of_array,dtype = mydtype)
        for j,thisname in enumerate(names):
            try:
                retval[thisname][:] = input[j][:]
            except:
                raise CustomError("error trying to load input for '"+thisname+"' of shape "+repr(shape(input[j]))+" into retval field of shape "+repr(shape(retval[thisname])))
        return retval
    else:
        try:
            return array([tuple(input)],dtype = mydtype)
        except:
            raise CustomError('problem trying to assign data of type',map(type,input),'\nvalues',input,'\nonto',mydtype,'\ndtype made from tuple:',zip(names,types,shapes))#,'"types" was',types)

#{{{ convert back and forth between lists, etc, and ndarray
def make_ndarray(array_to_conv,name_forprint = 'unknown',verbose = False): 
    if type(array_to_conv) in [int,int32,double,float,complex,complex128,float,bool,bool_]: # if it's a scalar
        pass
    elif type(array_to_conv) is str:
        pass
    elif type(array_to_conv) in [list,ndarray] and len(array_to_conv) > 0:
        array_to_conv = rec.fromarrays([array_to_conv],names = 'LISTELEMENTS') #list(rec.fromarrays([b])['f0']) to convert back
    elif type(array_to_conv) in [list,ndarray] and len(array_to_conv) is 0:
        array_to_conv = None
    elif array_to_conv is  None:
        pass
    else:
        raise CustomError('type of value (',type(array_to_conv),') for attribute name',name_forprint,'passed to make_ndarray is not currently supported')
    return array_to_conv

def unmake_ndarray(array_to_conv,name_forprint = 'unknown',verbose = False): 
    r'Convert this item to an ndarray'
    if (type(array_to_conv) is recarray) or (type(array_to_conv) is ndarray and array_to_conv.dtype.names is not None and len(array_to_conv.dtype.names)>0):
        #{{{ if it's a record/structured array, it should be either a list or dictionary
        if 'LISTELEMENTS' in array_to_conv.dtype.names:
            if array_to_conv.dtype.names == tuple(['LISTELEMENTS']):
                retval = list(array_to_conv['LISTELEMENTS'])
            else:
                raise CustomError('Attribute',name_forprint,'is a recordarray with a LISTELEMENTS field, but it also has other dimensions:',array_to_conv.dtype.names,'not',tuple(['LISTELEMENTS']))
        elif len(array_to_conv)==1:
            thisval = dict(zip(a.dtype.names,a.tolist()[0]))
        else: raise CustomError('You passed a structured array, but it has more than one dimension, which is not yet supported\nLater, this should be supported by returning a dictionary of arrays')
        #}}}
    elif type(array_to_conv) is ndarray and len(array_to_conv)==1:
        #{{{ if it's a length 1 ndarray, then return the element
        retval = array_to_conv.tolist()
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",type(array_to_conv),"is a numpy array of length one"
        #}}}
    elif type(array_to_conv) in [string_,int32,float64,bool_]:
        #{{{ map numpy strings onto normal strings
        retval = array_to_conv.tolist()
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",type(array_to_conv),"is a numpy scalar"
        #}}}
    elif type(array_to_conv) is list:
        #{{{ deal with lists
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"is a list"
        typeofall = map(type,array_to_conv)
        if all(map(lambda x: x is string_,typeofall)):
            if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",typeofall,"are all numpy strings"
            retval = map(str,array_to_conv)
        else:
            if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",typeofall,"are not all numpy string"
            retval = array_to_conv
        #}}}
    else:
        if verbose: print "(from unmake ndarray verbose):", name_forprint,"=",type(array_to_conv),"is not a numpy string or record array"
        retval = array_to_conv
    return retval

#}}}
#}}}

def emptytest(x): # test is it is one of various forms of empty
   if type(x) in [list,array]:
       if len(x) == 0:
           return True
       elif x is array(None):
           return True
       elif len(x) > 0:
           return False
       #don't want the following, because then I may need to pop, etc
       #if type(x) is list and all(map(lambda x: x is None,x)): return True
   if size(x) is 1 and x is None: return True
   if size(x) is 0: return True
   return False

def lsafen(*string,**kwargs):
    "see lsafe, but with an added double newline"
    string = list(string)
    string += ['\n\n']
    return lsafe(*tuple(string),**kwargs)

def lsafe(*string,**kwargs):
    "Output properly escaped for latex"
    if len(string) > 1:
        lsafewkargs = lambda x: lsafe(x,**kwargs)
        return ' '.join(list(map(lsafewkargs,string)))
    else:
        string = string[0]
    #{{{ kwargs
    spaces = False
    if 'spaces' in kwargs.keys():
        spaces = kwargs.pop('spaces')
    if 'wrap' in kwargs.keys():
        wrap = kwargs.pop('wrap')
    else:
        wrap = None
    #}}}
    if type(string) is not str:
        string = repr(string)
    if wrap is True:
        wrap = 60
    if wrap is not None:
        string = '\n'.join(textwrap.wrap(string,wrap))
    string = string.replace('\\','\\textbackslash ')
    if spaces:
        string = string.replace(' ','\\ ')
    string = string.replace('\n\t','\n\n\\quad ')
    string = string.replace('\t','\\quad ')
    string = string.replace('_',r'\_')
    string = string.replace('{',r'\{')
    string = string.replace('}',r'\}')
    string = string.replace('$$',r'ACTUALDOUBLEDOLLAR')
    string = string.replace(']',r'$]$')
    string = string.replace('[',r'$[$')
    string = string.replace('<',r'$<$')
    string = string.replace('>',r'$>$')
    string = string.replace('$$',r'')
    string = string.replace('ACTUALDOUBLEDOLLAR',r'$$')
    string = string.replace('^',r'\^')
    string = string.replace('#',r'\#')
    string = string.replace('%',r'\%')
    string = string.replace('&',r'\&')
    string = string.replace('+/-',r'\ensuremath{\pm}')
    string = string.replace('|',r'$|$')
    return string

#{{{ errors
class CustomError(Exception):
    def __init__(self, *value, **kwargs):
        if 'figurelist' in kwargs.keys():
            lplotfigures(kwargs.pop('figurelist'),'error_plots.pdf')
        if len(value)>1:
            retval = map(str,value)
        else:
            retval = str(value)
        retval = map(str,value)
        retval = ' '.join(retval)
        retval = '\n'+'\n'.join(textwrap.wrap(retval,90,replace_whitespace = False))
        if traceback.format_exc() != 'None':
            retval += '\n\nOriginal Traceback:\n'+''.join(['V']*40)+'\n'+traceback.format_exc() + '\n' + ''.join(['^']*40) + '\n'
        Exception.__init__(self,retval)
        return

def copy_maybe_none(input):
    if input == None:
        return None
    else:
        if type(input) is list:
            return map(copy,input)
        else:
            return input.copy()

def maprep(*mylist):
    mylist = list(mylist)
    for j in range(0,len(mylist)):
        if type(mylist[j]) is not str:
            mylist[j] = mylist[j].__repr__()
    return ' '.join(mylist)
#}}}

#{{{ HDF5 functions
#{{{ helper function for HDF5 search
def gensearch(labelname,format = '%0.3f',value = None,precision = None):
    'obsolete -- use h5gensearch'
    if value == None:
        raise CustomError('You must pass a value to gensearch')
    if precision == None:
        precision = value*0.01 # the precision is 1% of the value, if we don't give an argument
    searchstring_high = '(%s < %s + (%s))'%tuple([labelname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([labelname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return searchstring_low + ' & ' + searchstring_high
def h5searchstring(fieldname,value,format = '%g',precision = 0.01):
    'search AROUND a certain value (overcomes some type conversion issues) optional arguments are the format specifier and the fractional precision'
    precision *= value
    searchstring_high = '(%s < %s + (%s))'%tuple([fieldname]+[format]*2)
    #print "\n\nDEBUG check format:\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_high = searchstring_high%(value,precision)
    #print "\n\nDEBUG after substitution with",value,precision,":\\begin{verbatim}",searchstring_high,r'\end{verbatim}'
    searchstring_low = '(%s > %s - (%s))'%tuple([fieldname]+[format]*2)
    searchstring_low = searchstring_low%(value,precision)
    return '(' + searchstring_low + ' & ' + searchstring_high + ')'
#}}}
def h5loaddict(thisnode,verbose = False):
    #{{{ load all attributes of the node
    retval = dict([(x,thisnode._v_attrs.__getattribute__(x))
        for x in thisnode._v_attrs._f_list('user')])
    #}}}
    for k,v in retval.iteritems():#{{{ search for record arrays that represent normal lists
        retval[k]  = unmake_ndarray(v,name_forprint = k,verbose = verbose)
    if type(thisnode) is tables.table.Table:#{{{ load any table data
        if verbose: print "It's a table\n\n"
        if 'data' in retval.keys():
            raise CustomError('There\'s an attribute called data --> this should not happen!')
        retval.update({'data':thisnode.read()})
    elif type(thisnode) is tables.group.Group:
        #{{{ load any sub-nodes as dictionaries
        mychildren = thisnode._v_children
        for thischild in mychildren.keys():
            if thischild in retval.keys():
                raise CustomError('There\'s an attribute called ',thischild,' and also a sub-node called the',thischild,'--> this should not happen!')
            retval.update({thischild:h5loaddict(mychildren[thischild])})
        #}}}
    else:
        raise CustomError("I don't know what to do with this node:",thisnode)
    #}}}
    return retval
def h5child(thisnode,childname,clear = False,create = None,verbose = False):
    r'''grab the child, optionally clearing it and/or (by default) creating it'''
    #{{{ I can't create and clear at the same time
    if create and clear:
        raise CustomError("You can't call clear and create at the same time!\nJust call h5child twice, once with clear, once with create")
    if create is None:
        if clear == True:
            create = False
        else:
            create = True
    #}}}
    h5file = thisnode._v_file
    try:
        childnode = h5file.getNode(thisnode,childname)
        if verbose:
            print lsafe('found',childname)
        if clear:
            childnode._f_remove(recursive = True)
            childnode = None
    except tables.NoSuchNodeError:
        if create is False and not clear:
            raise CustomError('Trying to grab a node that does not exist with create = False')
        elif clear:
            childnode = None
        else:
            childnode = h5file.createGroup(thisnode,childname)
            if verbose:
                print lsafe('created',childname)
    return childnode
def h5remrows(bottomnode,tablename,searchstring):
    try:
        thistable = bottomnode.__getattr__(tablename)
        counter = 0
        try:
            data = thistable.readWhere(searchstring).copy()
        except:
            raise CustomError('Problem trying to remove rows using search string',searchstring,'in',thistable)
        for row in thistable.where(searchstring):
            if len(thistable) == 1:
                thistable.remove()
                counter += 1
            else:
                thistable.removeRows(row.nrow - counter,None) # counter accounts for rows I have already removed.
                counter += 1
        return counter,data
    except tables.NoSuchNodeError:
        return False,None
def h5addrow(bottomnode,tablename,*args,**kwargs):
    'add a row to a table, creating it if necessary, but don\'t add if the data matches the search condition'
    #{{{ process kwargs
    match_row = None
    if 'match_row' in kwargs.keys():
        match_row = kwargs.pop('match_row')
    verbose = False
    if 'verbose' in kwargs.keys():
        verbose = kwargs.pop('verbose')
    only_last = True
    if 'only_last' in kwargs.keys():
        only_last = kwargs.pop('only_last')
    if len(kwargs) != 0:
        raise ValueError('kwargs'+repr(kwargs)+'not understood!!')
    #}}}
    try: # see if the table exists
        mytable = h5table(bottomnode,tablename,None)
        #{{{ auto-increment "index"
        newindex = mytable.read()['index'].max() + 1L
        #}}}
        # here is where I would search for the existing data
        if match_row is not None:
            try:
                matches = mytable.readWhere(match_row)
            except NameError:
                raise CustomError('The columns available are',mytable.colnames)
            if len(matches) > 0:
                if only_last:
                    if verbose: print r'\o{',lsafen(len(matches),"rows match your search criterion, returning the last row"),'}'
                    return mytable,matches['index'][-1]
                else:
                    return mytable,matches['index'][:]
        tableexists = True
    except CustomError: # if table doesn't exist, create it
        newindex = 1L
        tableexists = False
    if len(args) == 1 and (type(args[0]) is dict):
        listofnames,listofdata = map(list,zip(*tuple(args[0].items())))
    elif len(args) == 2 and type(args[0]) is list and type(args[1]) is list:
        listofdata = args[0]
        listofnames = args[1]
    else:
        raise TypeError('h5addrow takes either a dictionary for the third argument or a list for the third and fourth arguments')
    try:
        listofdata = [newindex] + listofdata
    except:
        raise TypeError('newindex is'+repr(newindex)+'listofdata is'+repr(listofdata))
    listofnames = ['index'] + listofnames
    myrowdata = make_rec(listofdata,listofnames)
    if tableexists:
        try:
            mytable.append(myrowdata)
        except ValueError:
            print lsafen("I'm about to flag an error, but it looks like there was an issue appending",myrowdata)
            tabledforerr = mytable.read()
            #raise CustomError('Value of table',tabledforerr.dtype.,'while value of row',myrowdata.dtype)
            raise CustomError('Value of table -- compare names and values table data vs. the row you are trying to add\n','\n'.join(map(repr,zip(mytable.read().dtype.fields.keys(),
                        mytable.read().dtype.fields.values(),
                        myrowdata.dtype.fields.keys(),
                        myrowdata.dtype.fields.values()))))
            #raise CustomError('Value of table',mytable.read().dtype,'while value of row',myrowdata.dtype,'data in row:',myrowdata,
            #        'the one that doesn\'t match is',[mytable.read().dtype.descr[x]
            #            for x in list(myrowdata.dtype)
            #            if mytable.read().dtype.descr[x]==myrowdata.dtype.descr[x]])
        mytable.flush()
    else:
        recorddata = myrowdata
        try:
            mytable = h5table(bottomnode,
                    tablename,
                    recorddata)
        except:
            raise CustomError('Error trying to write record array:',repr(recorddata),'from listofdata',listofdata,'and names',listofnames)
        mytable.flush()
    return mytable,newindex
def h5table(bottomnode,tablename,tabledata):
    'create the table, or if tabledata is None, just check if it exists'
    #{{{ save but don't overwrite the table
    h5file = bottomnode._v_file
    if tablename not in bottomnode._v_children.keys():
        if tabledata is not None:
            datatable = h5file.createTable(bottomnode,tablename,tabledata) # actually write the data to the table
        else:
            raise CustomError('You passed no data, so I can\'t create table',tablename,'but it doesn\'t exist in',bottomnode,'which has children',bottomnode._v_children.keys())
    else:
        if tabledata is not None:
            raise CustomError('You\'re passing data to create the table, but the table already exists!')
        else:
            pass
    return bottomnode._v_children[tablename]
    #}}}
def h5nodebypath(h5path,verbose = False,force = False,only_lowest = False,check_only = False):
    r'''return the node based on an absolute path, including the filename'''
    if verbose: print lsafen("DEBUG: called h5nodebypath on",h5path)
    h5path = h5path.split('/')
    #{{{ open the file / check if it exists
    if verbose: print lsafen('h5path=',h5path)
    try:
        if h5path[0] in listdir('.'):
            if verbose: print 'DEBUG: file exists\n\n'
        else:
            if check_only: raise CustomError("You're checking for a node in a file that does not exist")
            if verbose: print 'DEBUG: file does not exist\n\n'
        mode = 'a'
        #if check_only: mode = 'r'
        h5file = tables.openFile(h5path[0],mode = mode,title = 'test file')
    except IOError:
        raise CustomError('I think the HDF5 file has not been created yet, and there is a bug pytables that makes it freak out, but you can just run again.')
    #}}}
    currentnode = h5file.getNode('/') # open the root node
    for pathlevel in range(1,len(h5path)):#{{{ step down the path
            clear = False
            create = True
            if only_lowest or check_only:
                create = False
            if pathlevel == len(h5path)-1: # the lowest level
                if only_lowest:
                    create = True
                if force:
                    clear = True
            safetoleaveopen = False
            try:
                currentnode = h5child(currentnode, # current node
                        h5path[pathlevel], # the child
                        verbose = verbose,
                        create = create,
                        clear = clear)
                if verbose: print lsafen("searching for node path: descended to node",currentnode)
            except:
                if verbose: print lsafen("searching for node path: got caught searching for node",h5path[pathlevel])
                h5file.close()
                #print lsafen("DEBUG: Yes, I closed the file")
                raise CustomError('Problem trying to load node ',h5path)
            #}}}
    return h5file,currentnode
def h5attachattributes(node,listofattributes,myvalues):
    #print "DEBUG 5: node passed to h5attachattributes",node
    if node is None:
        raise CustomError('Problem!, node passed to h5attachattributes: ',node,'is None!')
    h5file = node._v_file
    if isinstance(myvalues,nddata):
        attributevalues = map(lambda x: myvalues.__getattribute__(x),listofattributes)
    elif type(myvalues) is list:
        attributevalues = myvalues
    else:
        raise CustomError("I don't understand the type of myvalues, which much be a list or a nddata object, from which the attribute values are retrieved")
    listout = list(listofattributes)
    for j,thisattr in enumerate(listofattributes):
        thisval = attributevalues[j]
        if type(thisval) in [dict]:
            dictnode = h5child(node,
                    thisattr,
                    clear = True)
            dictnode = h5child(node,
                    thisattr,
                    create = True)
            h5attachattributes(dictnode,
                    thisval.keys(),
                    thisval.values())
            thisval = None
            listout.remove(thisattr)
        else:
            thisval = make_ndarray(thisval,name_forprint = thisattr)
        if thisval is not None:
            node._v_attrs.__setattr__(thisattr,thisval)
            listout.remove(thisattr)
    listofattributes[:] = listout # pointer
def h5inlist(columnname,mylist):
    'returns rows where the column named columnname is in the value of mylist'
    if type(mylist) is slice:
        if mylist.start is not None and mylist.stop is not None:
            return "(%s >= %g) & (%s < %g)"%(columnname,mylist.start,columnname,mylist.stop)
        elif mylist.stop is not None:
            return "(%s < %g)"%(columnname,mylist.stop)
        elif mylist.start is not None:
            return "(%s > %g)"%(columnname,mylist.start)
        else:
            raise ValueError()
    if type(mylist) is ndarray:
        mylist = mylist.tolist()
    if type(mylist) is not list:
        raise TypeError("the second argument to h5inlist must be a list!!!")
    if any([type(x) in [double,float64] for x in mylist]):
        if all([type(x) in [double,float64,int,int32,int64] for x in mylist]):
            return '('+'|'.join(map(lambda x: "(%s == %g)"%(columnname,x),mylist))+')'
    elif all([type(x) in [int,long,int32,int64] for x in mylist]):
        return '('+'|'.join(map(lambda x: "(%s == %g)"%(columnname,x),mylist))+')'
    elif all([type(x) is str for x in mylist]):
        return '('+'|'.join(map(lambda x: "(%s == '%s')"%(columnname,x),mylist))+')'
    else:
        raise TypeError("I can't figure out what to do with this list --> I know what to do with a list of numbers or a list of strings, but not a list of type"+repr(map(type,mylist)))
def h5join(firsttuple,secondtuple,
    additional_search = '',
    select_fields = None,
    pop_fields = None,
    verbose = False):
    #{{{ process the first argument as the hdf5 table and indices, and process the second one as the structured array to join onto
    if not ((type(firsttuple) is tuple) and (type(secondtuple) is tuple)):
        raise ValueError('both the first and second arguments must be tuples!')
    if not ((len(firsttuple) == 2) and (len(secondtuple) == 2)):
        raise ValueError('The length of the first and second arguments must be two!')
    tablenode = firsttuple[0]
    tableindices = firsttuple[1]
    if verbose: print 'h5join tableindices looks like this:',tableindices
    if type(tableindices) is not list:
        tableindices = [tableindices]
    if verbose: print 'h5join tableindices looks like this:',tableindices
    mystructarray = secondtuple[0].copy()
    mystructarrayindices = secondtuple[1]
    if type(mystructarrayindices) is not list:
        mystructarrayindices = [mystructarrayindices]
    #}}}
    #{{{ generate a search string to match potentially more than one key
    search_string = []
    if len(tableindices) != len(mystructarrayindices):
        raise ValueError('You must pass either a string or a list for the second element of each tuple!\nIf you pass a list, they must be of the same length, since the field names need to line up!')
    # this can't use h5inlist, because the and needs to be on the inside
    #{{{ this loop creates a list of lists, where the inner lists are a set of conditions that need to be satisfied
    # this is actually not causing  any trouble right now, but needs to be fixed, because of the way that it's doing the type conversion
    for thistableindex,thisstructarrayindex in zip(tableindices,mystructarrayindices):
        if thisstructarrayindex not in mystructarray.dtype.names:
            raise ValueError(repr(thisstructarrayindex)+" is not in "+repr(mystructarray.dtype.names))
        if type(mystructarray[thisstructarrayindex][0]) in [str,str_]:
            search_string.append(["(%s == '%s')"%(thistableindex,x) for x in mystructarray[thisstructarrayindex]])
        elif type(mystructarray[thisstructarrayindex][0]) in [int,double,float,float64,float32,int32,int64]:
            search_string.append(["(%s == %s)"%(thistableindex,str(x)) for x in mystructarray[thisstructarrayindex]])
            #print 'a g mapping for',[x for x in mystructarray[thisstructarrayindex]],'gives',search_string[-1],'\n\n'
        else:
            raise TypeError("I don't know what to do with a structured array that has a row of type"+repr(type(mystructarray[thisstructarrayindex][0])))
    #}}}
    search_string = [' & '.join(x) for x in zip(*tuple(search_string))] # this "and"s together the inner lists, since all conditions must be matched
    search_string = '('+'|'.join(search_string)+')' # then, it "or"s the outer lists, since I want to collect data for all rows of the table
    #}}}
    if len(additional_search) > 0:
        additional_search = " & (%s)"%additional_search
        search_string = search_string + additional_search
    if verbose: print '\n\nh5join generated the search string:',lsafen(search_string)
    retval = tablenode.readWhere(search_string)
    #{{{ then join the data together
    # here I'm debugging the join function, again, and again, and agin
    try:
        retval = decorate_rec((retval,tableindices),(mystructarray,mystructarrayindices)) # this must be the problem, since the above looks fine
    except:
        raise CustomError('Some problems trying to decorate the table',retval,'of dtype',retval.dtype,'with the structured array',mystructarray,'of dtype',mystructarray.dtype)
    if pop_fields is not None:
        if select_fields is not None:
            raise ValueError("It doesn't make sense to specify pop_fields and select_fields at the same time!!")
        select_fields = list(set(retval.dtype.names) ^ set(pop_fields))
    if select_fields is not None:
        if verbose: print '\n\nh5join original indices',lsafen(retval.dtype.names)
        try:
            retval = retval[select_fields]
        except ValueError:
            raise CustomError('One of the fields',select_fields,'is not in',retval.dtype.names)
    #}}}
    return retval
#}}}
#{{{ indices to slice
#}}}
#{{{ add slashes for dir's
def dirformat(file):
        #{{{ format strings
        if file[-1]!='/':
            file += '/'
        #}}}
        return file
#}}}
#}}}




#{{{general functions
def box_muller(length):
    r'''algorithm to generate normally distributed noise'''
    s1 = rand(length)
    s2 = rand(length)
    n1 = sqrt(-2*log(s1))*cos(2*pi*s2)
    n2 = sqrt(-2*log(s1))*sin(2*pi*s2)
    return (n1 + 1j * n2)*0.5
#}}}



def fa(input,dtype='complex128'):# make a fortran array
    return array(input,order='F',dtype=dtype) # will need transpose reverses the dimensions, since the bracketing still works in C order (inner is last index), but F tells it to store it appropriately in memory

def ndgrid(*input):
    thissize = list([1])
    thissize = thissize * len(input)
    output = list()
    for j in range(0,len(input)):
        tempsize = copy(thissize)
        tempsize[j] = input[j].size
        output.append(input[j].reshape(tempsize))
    return output
def pinvr(C,alpha):
    U,S,V = svd(C,full_matrices=0)
    #print 'U S V shapes:'
    #print U.shape
    #print S.shape
    #print V.shape
    if any(~isfinite(U)):
        raise CustomError('pinvr error, U is not finite')
    if any(~isfinite(V)):
        raise CustomError('pinvr error, V is not finite')
    if any(~isfinite(S)):
        raise CustomError('pinvr error, S is not finite')
    S = diag(S / (S**2 + alpha**2))
    if any(~isfinite(S)):
        raise CustomError('pinvr error, problem with S/(S^2+alpha^2) --> set your regularization higher')
    return dot(conj(transpose(V)),
            dot(S,conj(transpose(U))))
def sech(x):
    return 1./cosh(x)

def myfilter(x,center = 250e3,sigma = 100e3): ###???
    x = (x-center)**2
    x /= sigma**2
    return exp(-x)
#}}}

def sqrt(arg):
    if isinstance(arg,nddata):
        return arg**0.5
    elif isinstance(arg,sympy.symbol.Symbol):
        return sympy.sqrt(arg)
    else:
        return np_sqrt(arg)

if paramset.myparams['figlist_type'] == 'figlistl':
    from fornotebook import *
    figlist_var = figlistl
elif paramset.myparams['figlist_type'] == 'figlist':
    def obsn(*x): #because this is used in fornotebook, and I want it defined
        print ''.join(x),'\n'
    def obs(*x): #because this is used in fornotebook, and I want it defined
        print ''.join(map(repr,x))
    def lrecordarray(*x,**kwargs):
        print x # if I'm not using tex, it's easier to not use the formatting
    figlist_var = figlist
    def lsafe(*string,**kwargs):
        "replacement for normal lsafe -- no escaping"
        if len(string) > 1:
            lsafewkargs = lambda x: lsafe(x,**kwargs)
            return ' '.join(list(map(lsafewkargs,string)))
        else:
            string = string[0]
        #{{{ kwargs
        spaces = False
        if 'spaces' in kwargs.keys():
            spaces = kwargs.pop('spaces')
        if 'wrap' in kwargs.keys():
            wrap = kwargs.pop('wrap')
        else:
            wrap = None
        #}}}
        if type(string) is not str:
            string = repr(string)
        if wrap is True:
            wrap = 60
        if wrap is not None:
            string = '\n'.join(textwrap.wrap(string,wrap))
        return string
