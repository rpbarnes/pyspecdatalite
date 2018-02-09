
#{{{nddata
#{{{ shaping and allocating
class ndshape ():
    def __init__(self,*args):
        self.zero_dimensional = False
        if len(args) == 2:
            self.shape = list(args[0])
            self.dimlabels = args[1]
        if len(args) == 1: #assum that it's an nddata object
            if isinstance(args[0],nddata):
                self.shape = list(args[0].data.shape)
                self.dimlabels = list(args[0].dimlabels)
                if len(self.shape) == 0 and len(self.dimlabels) == 0:
                    self.zero_dimensional = True
                    return
            else:
                raise ValueError('If you pass a single argument, it must be an nddata')
        return

    def __setitem__(self,reference,setto):
        self.shape[self.axn(reference)] = setto
        return

    def axn(self,axis):
        r'return the number for the axis with the name "axis"'
        try:
            return self.dimlabels.index(axis)
        except:
            raise CustomError('there is no axis named',axis,'all axes are named',self.dimlabels)

    def copy(self):
        try:
            return deepcopy(self)
        except:
            raise RuntimeError('Some type of error trying to run deepcopy on'+repr(self))

    def matchdims(self,arg):
        r'returns shape with [not in self, len 1] + [overlapping dims between arg + self] + [not in arg] --> this is better accomplished by using sets as I do in the matchdims below'
        for k in set(self.dimlabels) & set(arg.dimlabels):
            a = arg.shape[arg.axn(k)]
            b = self.shape[self.axn(k)]
            if a != b:
                raise CustomError('the',k,'dimension is not the same for self',self,'and arg',arg)
        if isinstance(arg,nddata):
            arg = ndshape(arg)
        #{{{ add extra 1-len dims
        addeddims = set(self.dimlabels) ^ set(arg.dimlabels) & set(arg.dimlabels)
        self.dimlabels = list(addeddims) + self.dimlabels
        self.shape = [1] * len(addeddims) + list(self.shape)
        #}}}
        return self

    def __add__(self,arg):
        'take list of shape,dimlabels'
        shape = arg[0]
        dimlabels = arg[1]
        if type(shape) is str:
            shape,dimlabels = dimlabels,shape
        if isscalar(self.shape):
            self.shape = [self.shape]
        if isscalar(self.dimlabels):
            self.dimlabels = [self.dimlabels]
        if isscalar(shape):
            shape = [shape]
        if isscalar(dimlabels):
            dimlabels = [dimlabels]
        self.shape = shape + self.shape
        self.dimlabels = dimlabels + self.dimlabels
        return self

    def add_correctly(self,arg):
        '''take list of shape,dimlabels
        this is the correct function, until I can fix my back-references for add, which does it backwards'''
        shape = arg[0]
        dimlabels = arg[1]
        if type(shape) is str:
            shape,dimlabels = dimlabels,shape
        if isscalar(self.shape):
            self.shape = [self.shape]
        if isscalar(self.dimlabels):
            self.dimlabels = [self.dimlabels]
        if isscalar(shape):
            shape = [shape]
        if isscalar(dimlabels):
            dimlabels = [dimlabels]
        self.shape = self.shape + shape
        self.dimlabels = self.dimlabels + dimlabels
        return self

    def __repr__(self): #how it responds to print
        return zip(self.shape,self.dimlabels).__repr__()
    
    def __getitem__(self,args):
        try:
            mydict = dict(zip(self.dimlabels,self.shape))
        except:
            raise CustomError("either dimlabels=",self.dimlabels,"or shape",self.shape,"not in the correct format")
        try:
            return mydict[args]
        except:
            raise CustomError("one or more of the dimensions named",args,"do not exist in",self.dimlabels)

    def pop(self,label):
        r'remove a dimension'
        thisindex = self.axn(label)
        self.dimlabels.pop(thisindex)
        self.shape.pop(thisindex)
        return self

    def alloc(self,dtype='complex128',labels = False,format = 0):
        try:
            if format == 0:
                try:
                    emptyar = zeros(tuple(self.shape),dtype=dtype)
                except TypeError:
                    raise TypeError("You passed a type of "+repr(dtype)+", which was likely not understood (you also passed a shape of "+repr(tuple(self.shape))+")")
            elif format == 1:
                emptyar = ones(tuple(self.shape),dtype=dtype)
            elif format is None:
                emptyar = empty(tuple(self.shape),dtype=dtype)
            else:
                emptyar = format*ones(tuple(self.shape),dtype=dtype)
        except TypeError:
            raise CustomError('Wrong type for self.shape',map(type,self.shape))
        retval = nddata(emptyar,self.shape,self.dimlabels)
        if labels:
            retval.labels(self.dimlabels,map(lambda x: double(r_[0:x]),self.shape))
        return retval
#}}}
#{{{ format out to a certain decimal place
def dp(number,decimalplaces,scientific=False):
    if scientific:
        tenlog = floor(log(number)/log(10.))
        number /= 10**tenlog
        fstring = '%0.'+'%d'%decimalplaces+r'f\times 10^{%d}'%tenlog
    else:
        fstring = '%0.'+'%d'%decimalplaces+'f'
    return fstring%number
#}}}
#{{{ concatenate datalist along dimname
def concat(datalist,dimname,chop = False,verbose = False):
    #{{{ allocate a new datalist structure  
    newdimsize = 0
    #print 'DEBUG: type(datalist)',type(datalist)
    try:
        shapes = map(ndshape,datalist)
    except:
        if type(datalist) is not list:
            raise CustomError('You didn\'t pass a list, you passed a',type(datalist))
        raise CustomError('Problem with what you passed to concat, list of types,',map(type,datalist))
    other_info_out = datalist[0].other_info
    for j in range(0,len(datalist)):
        #{{{ make list for the shape to check, which contains the dimensions we are NOT concatting along
        if dimname in shapes[j].dimlabels:
            newdimsize += shapes[j][dimname]
            shapetocheck = list(shapes[j].shape)
            shapetocheck.pop(shapes[j].axn(dimname))
        else:
            newdimsize += 1
            shapetocheck = list(shapes[j].shape)
        #}}}
        if j is 0:
            shapetocheckagainst = shapetocheck
        else:
            if any(~(array(shapetocheck) == array(shapetocheckagainst))):
                if chop:
                    if verbose:
                        print lsafen(repr(shapetocheck)),lsafen(repr(shapetocheckagainst))
                        raise CustomError('For item ',j,'in concat, ',shapetocheck,'!=',shapetocheckagainst,'where all the shapes of the things you\'re trying to concat are:',shapes)
                else:
                    raise CustomError('For item ',j,'in concat, ',shapetocheck,'!=',shapetocheckagainst,'where all the shapes of the things you\'re trying to concat are:',shapes)
    newdatalist = ndshape(datalist[-1])
    if dimname in newdatalist.dimlabels:
        newdatalist[dimname] = newdimsize
    else:
        newdatalist += ([newdimsize],[dimname])
    #print "DEBUG newdatalist is shaped like",newdatalist
    try:
        newdatalist = newdatalist.alloc()
    except:
        raise CustomError("trying to alloc the newdatalist",newdatalist,"created a problem")
    if datalist[0].get_error() is not None:
        newdatalist.set_error(zeros(shape(newdatalist.data)))
    #}}}
    #{{{ actually contract the datalist
    newdimsize = 0 # now use it to track to position
    for j in range(0,len(datalist)):
        if dimname in shapes[j].dimlabels:
            newdatalist[dimname,newdimsize:newdimsize+shapes[j][dimname]] = datalist[j]
            newdimsize += shapes[j][dimname]
        else:
            newdatalist[dimname,newdimsize:newdimsize+1] = datalist[j]
            newdimsize += 1
    #}}}
    #{{{ pull the axis labels from the last item in the list
    if len(datalist[-1].axis_coords)>0:
        dimlabels = list(datalist[-1].dimlabels)
        axis_coords = list(datalist[-1].axis_coords)
        #print "axis_coords are",axis_coords,"for",dimlabels
        if dimname in dimlabels:
            thisindex = dimlabels.index(dimname)
            dimlabels.pop(thisindex)
            axis_coords.pop(thisindex)
        dimlabels += [dimname]
        axis_coords += [r_[0:newdimsize]]
        try:
            newdatalist.labels(dimlabels,axis_coords)
        except:
            raise CustomError("trying to attach axes of lengths",map(len,axis_coords),"to",dimlabels)
    #}}}
    newdatalist.other_info = other_info_out
    return newdatalist
#}}}

class nddata (object):
    want_to_prospa_decim_correct = False
    def __init__(self,*args,**kwargs):
        if len(args) > 1:
            if len(args) == 2:
                if len(args[0].shape) == 1 and type(args[1]) is str:
                    self.__my_init__(args[0],[len(args[0])],[args[1]])
                    self.labels(args[1],args[0])
                else:
                    raise ValueError('You can pass two arguments only if you pass a 1d ndarray and a name for the axis') 
            else:
                self.__my_init__(args[0],args[1],args[2],**kwargs)
        else:
            self.__my_init__(args[0],[-1],['value'],**kwargs)
        return

    def __my_init__(self,data,sizes,dimlabels,axis_coords=[],ft_start_time = 0.,data_error = None, axis_coords_error = None,axis_coords_units = None, data_units = None, other_info = {}):
        self.genftpairs = False
        if not (type(data) is ndarray):
            #if (type(data) is float64) or (type(data) is complex128) or (type(data) is list):
            if isscalar(data) or (type(data) is list) or (type(data) is tuple):
                data = array(data)
            else:
                raise CustomError('data is not an array, it\'s',type(data),'!')
        if not (type(dimlabels) is list):
            raise CustomError('labels are not a list')
        try:
            self.data = reshape(data,sizes)
        except:
            raise CustomError("trying to reshape a ",data.shape,"array with list of sizes",sizes)
        self.dimlabels = dimlabels
        self.axis_coords = axis_coords
        #if len(axis_coords) > 0:
        #    testshape = data.shape
        #    if not all([len(axis_coords[j])==testshape[j] if axis_coords[j] is not None else True for j in range(0,len(axis_coords))]):
        #        raise IndexError('The length of your axis labels (axis_coords) (shape %s) and your axis data (shape %s) does not match!!!'%(repr([len(thiscoord) for thiscoord in axis_coords]),repr(data.shape)))
        self.ft_start_time = ft_start_time
        self.data_error = data_error
        self.data_units = data_units
        self.other_info = dict(other_info)
        if axis_coords_error == None:
            self.axis_coords_error = [None]*len(axis_coords)
        else:
            self.axis_coords_error = axis_coords_error
        if axis_coords_units == None:
            self.axis_coords_units = [None]*len(axis_coords)
        else:
            self.axis_coords_units = axis_coords_units 
        return

    def _contains_symbolic(self,string):
        return string[:9] == 'symbolic_' and hasattr(self,string)
    #{{{ for printing
    def __repr__(self):
        retval = repr(self.data) 
        retval += '\n\t\t+/-'
        retval += repr(self.get_error())
        retval += '\n\tdimlabels=['
        retval += repr(self.dimlabels)
        retval += ']\n\taxes='
        def rep_this_dict(starting_indent,thisdict,errordict):
            dictrep = []
            for k,v in thisdict.iteritems():
                dictrep.append('`'+k+'\':'+repr(v)+starting_indent+'\t\t+/-'+repr(errordict[k]))
            return '{'+(','+starting_indent+'\t').join(dictrep)+'}' # separate with an extra comma, the existing indent, and a tab
        retval += rep_this_dict('\n\t',self.mkd(self.axis_coords),self.mkd(self.axis_coords_error))
        #retval += '\n\t\t+/-'
        #retval += rep_this_dict('\n\t\t',self.mkd(self.axis_coords_error))
        retval += '\n'
        return retval
    #}}}
    #{{{ for plotting

    def gnuplot_save(self,filename):
        x = self.getaxis(self.dimlabels[0])[:5]
        y = self.getaxis(self.dimlabels[1])[:5]
        z = self.data[:5,:5]
        print "size of x",size(x),"size of y",size(y),"size of z",size(z)
        print "x",x,"y",y,"z",z
        data = empty((z.shape[0]+1,z.shape[1]+1))
        data[1:,1:] = z[:]
        data[0,0] = z.shape[1]
        data[0,1:] = y.flatten()
        data[1:,0] = x.flatten()
        print "data",data
        fp = open('auto_figures/'+filename+'.dat','w')
        fp.write(float32(data).tostring())
        fp.write('\n')
        fp.close()
        return
    #{{{ sort and shape the data for 3d plotting
    def sorted_and_xy(self):
        sortedself = self.copy()
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if len(sortedself.dimlabels) > 2:
            raise CustomError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = sortedself.dimlabels[0]
        y_dim = sortedself.dimlabels[1]
        x_axis = sortedself.retaxis(x_dim).data
        y_axis = sortedself.retaxis(y_dim).data
        #}}}
        return sortedself,x_axis,y_axis
    def matrices_3d(self,also1d = False,invert = False):
        sortedself,x_axis,y_axis = self.sorted_and_xy()
        if invert:
            print "trying to invert meshplot-like data"
        X = x_axis*ones(shape(y_axis))
        Y = ones(shape(x_axis))*y_axis
        Z = real(sortedself.data)
        if invert:
            X = X[:,::-1]
            Y = Y[:,::-1]
            Z = Z[:,::-1]
        if also1d:
            if invert:
                return X,Y,Z,x_axis[::-1],y_axis[::-1]
            else:
                return X,Y,Z,x_axis,y_axis
        else:
            return X,Y,Z
    #}}}

    def mayavi_surf(self):
        from mayavi.mlab import surf
        X,Y,Z = self.matrices_3d()
        s = surf(X,Y,Z)
        return s
    #}}}

    #{{{ error-related functions
    def normalize(self,axis,first_figure = None):#,whichpoint = slice(0,1,None)):
       x = self.data
       n = len(x)
       S = sparse.lil_matrix((n,n))
       S.setdiag((self.get_error())**2)
       self.set_error(None)
       first_point = self[axis,0:1].copy() # this makes another instance that contains just the first point, for error propagation
       B = sparse.lil_matrix((n,n))
       B.setdiag(1./x)
       B[0,:] = -x/(x[0]**2) # Sparse seems to only support row assignment, so make the transpose to give it what it wants
       B[0,0] = 0.0
       B = B.T
       E = B * S * (B.T) # verified that this is matrix multiplication
       self /= first_point # this gives the experimentally measured E
       #{{{ now, chop out the first point, which is meaningless
       self = self[axis,1:]
       E = E[1:,1:]
       #}}}
       self.set_error(sqrt(E.diagonal()))
       E.setdiag(zeros(n-1))
       self.data_covariance = E
       return self

    def get_covariance(self):
        '''this returns the covariance matrix of the data'''
        if hasattr(self,'data_covariance'):
            E = self.data_covariance.copy()
        else:
            n = size(self.data)
            E = sparse.lil_matrix((n,n))
        try:
            E.setdiag(self.get_error()**2)
        except:
            raise CustomError('Problem getting error because error is',self.get_error())
        return E.toarray()
    #}}}

    #{{{ shortcuts for axes
    def axlen(self,axis):
        return shape(self.data)[self.axn(axis)]
    def axn(self,axis):
        r'return the number for the axis with the name "axis"'
        try:
            return self.dimlabels.index(axis)
        except:
            raise CustomError('there is no axis named',axis,'all axes are named',self.dimlabels)
    #}}}

    #{{{ dictionary functions -- these convert between two formats:
    # dictionary -- stuff labeled according the dimension label.
    # list -- same information, but it's assumed they are listed in the order given by "dimlabels"
    def mkd(self,*arg,**kwargs):
        'make dictionary format'
        #{{{ process kwargs
        give_None = True
        if len(kwargs) > 0:
            if 'give_None' in kwargs.keys():
                give_None = kwargs.pop('give_None')
        if len(kwargs) > 0:
            raise CustomError("you passed mkd kwargs I didn't understand:",kwargs)
        #}}}
        if len(arg) == 1:
            if emptytest(arg[0]):
                return dict(zip(self.dimlabels,
                    [None]*len(self.dimlabels)))
            if len(arg[0]) != len(self.dimlabels):
                print r"{\color{red}WARNING! mkd error (John will fix this later):}"
                print "When making a dictionary with mkd, you must pass a list that has one element for each dimension!  dimlabels is "+repr(self.dimlabels)+" and you passed "+repr(arg)+'\n\n'
                raise ValueError("When making a dictionary with mkd, you must pass a list that has one element for each dimension!  dimlabels is "+repr(self.dimlabels)+" and you passed "+repr(arg))
            for i,v in enumerate(arg[0]):
                if type(v) == ndarray:
                    if v.shape == ():
                        arg[0][i] = None
            if give_None:
                return dict(zip(self.dimlabels,arg[0]))
            else:
                #{{{ don't return values for the things that are None
                mykeys = [self.dimlabels[j] for j in range(0,len(self.dimlabels)) if arg[0][j] is not None]
                myvals = [arg[0][j] for j in range(0,len(self.dimlabels)) if arg[0][j] is not None]
                return dict(zip(mykeys,myvals))
                #}}}
        elif len(arg) == 0:
            if not give_None:
                raise CustomError("You can't tell me not to give none and then not pass me anything!!")
            return dict(zip(self.dimlabels,
                [None]*len(self.dimlabels)))
        else:
            raise CustomError('.mkd() doesn\'t know what to do with %d arguments',len(arg))

    def fld(self,dict_in,noscalar = False):
        'flatten dictionary -- return list'
        return [dict_in[x] for x in self.dimlabels]
    #}}}

    #{{{ set + get the error + units
    #{{{ set units
    def set_units(self,*args):
        if len(args) == 2:
            unitval = args[1] # later, have some type of processing bojive
            if self.axis_coords_units == None or len(self.axis_coords_units) == 0:
                self.axis_coords_units = [None] * len(self.dimlabels)
            self.axis_coords_units[self.axn(args[0])] = unitval
        elif len(args) == 1:
            unitval = args[0] # later, have some type of processing bojive
            self.data_units = unitval
        else:
            raise CustomError(".set_units() takes data units or 'axis' and axis units")
        return self

    def human_units(self,verbose = False):
        prev_label = self.get_units()
        oom_names =   ['T' , 'G' , 'M' , 'k' , '' , 'm' , '\\mu ' , 'n' , 'p']
        oom_values = r_[12 , 9   , 6   , 3   , 0  , -3  , -6     , -9  , -12]
        if prev_label is not None and len(prev_label)>0:
            #{{{ find the average order of magnitude, rounded down to the nearest power of 3
            average_oom = log10(abs(self.data))/3.
            average_oom = average_oom[isfinite(average_oom)].mean()
            #}}}
            if verbose: print "(human units): for data the average oom is",average_oom*3
            if round(average_oom) == 0.0:
                average_oom = 0
            else:
                average_oom = 3*floor(average_oom)
            if verbose: print "(human units): for data I round this to",average_oom
            oom_index = argmin(abs(average_oom-oom_values))
            if verbose: print "(human units): for data, selected an oom value of",oom_values[oom_index]
            self.data[:] /= 10**oom_values[oom_index]
            self.set_units(oom_names[oom_index]+prev_label)
        else:
            if verbose: print 'data does not have a unit label'
        for thisaxis in self.dimlabels:
            prev_label = self.get_units(thisaxis)
            if prev_label is not None and len(prev_label)>0:
                data_to_test = self.getaxis(thisaxis)
                if data_to_test is not None:
                    try:
                        data_to_test = data_to_test[isfinite(data_to_test)]
                    except:
                        raise CustomError('data_to_test is',data_to_test,'isfinite is',isfinite(data_to_test))
                    #{{{ find the average order of magnitude, rounded down to the nearest power of 3
                    average_oom = log10(abs(data_to_test))/3.
                    average_oom = average_oom[isfinite(average_oom)].mean()
                    #}}}
                    if verbose: print "(human units): for axis",thisaxis,"the average oom is",average_oom*3
                    average_oom = 3*floor(average_oom)
                    if verbose: print "(human units): for axis",thisaxis,"I round this to",average_oom
                    oom_index = argmin(abs(average_oom-oom_values))
                    if verbose: print "(human units): for axis",thisaxis,"selected an oom value of",oom_values[oom_index]
                    x = self.getaxis(thisaxis)
                    x[:] /= 10**oom_values[oom_index]
                    self.setaxis(thisaxis,x)
                    self.set_units(thisaxis,oom_names[oom_index]+prev_label)
                else:
                    if verbose: print thisaxis,'does not have an axis label'
            else:
                if verbose: print thisaxis,'does not have a unit label'
        return self
    #}}}

    #{{{ get units
    def units_texsafe(self,*args):
        retval = self.get_units(*args)
        if retval is None:
            return None
        if retval.find('\\') > -1:
            retval = '$'+retval+'$'
        return retval

    def replicate_units(self,other):
        for thisaxis in self.dimlabels:
            if other.get_units(thisaxis) is not None:
                self.set_units(thisaxis,other.get_units(thisaxis))
        if other.get_units() is not None:
            self.set_units(other.get_units(thisaxis))
        return self

    def get_units(self,*args):
        if len(args) == 1:
            if self.axis_coords_units == None:
                return None
            if len(self.axis_coords_units) == 0:
                return None
            try:
                return self.axis_coords_units[self.axn(args[0])]
            except:
                raise CustomError('problem getting units for',args[0],'dimension',self.dimlabels,self.axis_coords_units)
        elif len(args) == 0:
            return self.data_units
        else:
            raise CustomError(".set_units() takes axis or nothing")
    #}}}

    #{{{ set error
    def set_error(self,*args):
        '''set the errors\neither set_error('axisname',error_for_axis) or set_error(error_for_data)'''
        if (len(args) is 1) and isscalar(args[0]):
            args = (r_[args[0]],)
        if (len(args) is 1) and (type(args[0]) is ndarray):
            self.data_error = reshape(args[0],shape(self.data))
        elif (len(args) is 1) and (type(args[0]) is list):
            self.data_error = reshape(array(args[0]),shape(self.data))
        elif (len(args) is 2) and (type(args[0]) is str) and (type(args[1]) is ndarray):
            self.axis_coords_error[self.axn(args[0])] = args[1]
        elif (len(args) is 2) and (type(args[0]) is str) and (args[1] is None):
            self.axis_coords_error[self.axn(args[0])] = args[1]
        elif (len(args) is 1) and args[0] is None:
            self.data_error = None
        else:
            raise CustomError('Not a valid argument to set_error:',map(type,args))
        return self
    #}}}

    #{{{ random mask -- throw out points
    def random_mask(self,axisname,threshold = exp(-1.0),inversion = False):
        r'''generate a random mask with about 'threshold' of the points thrown out'''
        if inversion:
            threshold = threshold / (1.0 - threshold)
        myr = rand(self.data.shape[self.axn(axisname)]) # random array same length as the axis
        return myr > threshold
    #}}}

    #{{{ get error
    def get_error(self,*args):
        '''get a copy of the errors\neither set_error('axisname',error_for_axis) or set_error(error_for_data)'''
        if (len(args) is 0):
            if self.data_error is None:
                return None
            else:
                return real(self.data_error)
        elif (len(args) is 1):
            thearg = args[0]
            if type(thearg) is str_:
                thearg = str(thearg) # like in the other spot, this became necessary with some upgrade, though I'm not sure that I should maybe just change the error functions to treat the numpy string in the same way
            if (type(thearg) is str):
                if len(self.axis_coords_error) == 0: self.axis_coords_error = [None] * len(self.dimlabels) # is we have an empty axis_coords_error, need to fill with None's
                try:
                    errorforthisaxis = self.axis_coords_error[self.axn(thearg)]
                except:
                    raise CustomError('Problem trying to load error',self.axn(thearg),'for axis',thearg,'out of',self.axis_coords_error)
                if errorforthisaxis is None:
                    return None
                else:
                    x = self.axis_coords_error[self.axn(thearg)]
                    if type(x) is ndarray:
                        if x.shape == ():
                            return None
                        else:
                            return real(self.axis_coords_error[self.axn(thearg)])
                    else:
                        return real(self.axis_coords_error[self.axn(thearg)])
        else:
            raise CustomError('Not a valid argument to get_error: *args=',args,'map(type,args)=',map(type,args))
        #}}}
    #}}}

    #{{{ match dims --
    def matchdims(self,other):
        r'add any dimensions to self that are not present in other'
        #print 'diagnose: matching',ndshape(self),'to',ndshape(other)
        addeddims =  list(set(self.dimlabels)^set(other.dimlabels))
        newdims = addeddims + self.dimlabels
        newshape = [1]*len(addeddims) + list(self.data.shape)
        #print 'diagnose: newshape',newshape,'newdims',newdims
        #{{{ reshape to the new dimensions  
        new_axis_coords = [r_[1]]*len(addeddims) + self.axis_coords
        self.data = self.data.reshape(newshape)
        self.dimlabels = newdims
        if len(self.axis_coords)>0:
            self.axis_coords = new_axis_coords
        #}}}
        #{{{ if we are adding dimensions, we will need to reorder to match the order of the other   
        if len(addeddims)>0:
            self.reorder(other.dimlabels)
        #}}}
        return self
    #}}}

    #{{{ rename
    def rename(self,previous,new):
        self.dimlabels[self.dimlabels.index(previous)] = new
        return self
    #}}}
    #{{{ display and other properties
    #{{{ set and get prop
    def set_prop(self,propname,val):
        self.other_info.update({propname:val})
        return
    def get_prop(self,propname):
        if propname not in self.other_info.keys():
            return None
        return self.other_info[propname]
    def name(self,*arg):
        r"""args:
           .name(newname) --> Name the object (for storage, etc)
           .name() --> Return the name"""
        if len(arg) == 1:
            self.set_prop('name',arg[0])
        elif len(arg) == 0:
            return self.get_prop('name')
        else:
            raise ValueError("invalid number of arguments")
    #}}}
    #{{{ set and get plot color
    def set_plot_color(self,thiscolor):
        if thiscolor is None:
            return
        if thiscolor is str:
            colordict = {'r':[1,0,0],
                    'g':[0,1,0],
                    'b':[0,0,1],
                    'k':[0,0,0],
                    'y':[0.5,0.5,0],
                    'o':[0.75,0.25,0],
                    'c':[0,0.5,0.5]}
            try:
                thiscolor = colordict[thiscolor]
            except:
                raise CustomError('Color',thiscolor,'not in dictionary')
        self.other_info.update({'plot_color':thiscolor})
        return
    def get_plot_color(self):
        if 'plot_color' in self.other_info.keys():
            return self.other_info['plot_color']
        else:
            return None
    #}}}
    #}}}
    #{{{ arithmetic
    def __add__(self,arg):
        if isscalar(arg):
            A = self.copy()
            if type(arg) is complex and self.data.dtype not in [complex128,complex64]:
                A.data = complex128(A.data)
            A.data += arg
            # error does not change
            return A
        #{{{ shape and add
        A,B = self.aligndata(arg)
        retval = A.copy()
        retval.data += B.data
        #}}}
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0
        if Aerr != None:
            Rerr += (Aerr)**2
        if Berr != None:
            Rerr += (Berr)**2
        Rerr = sqrt(real(Rerr)) # convert back to stdev
        if Aerr == None and Berr == None:
            Rerr = None
        retval.set_error(Rerr)
        return retval
    def __sub__(self,arg):
        return self.__add__(-1*arg)
    def __lt__(self,arg):
        if type(arg) is ndarray:
            retval = self.copy()
            retval.data = retval.data < arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data < B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data < arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __gt__(self,arg):
        if type(arg) is ndarray:
            retval = self.copy()
            retval.data = retval.data > arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data > B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data > arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __le__(self,arg):
        if type(arg) is ndarray:
            retval = self.copy()
            retval.data = retval.data <= arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data <= B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data <= arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __ge__(self,arg):
        if type(arg) is ndarray:
            retval = self.copy()
            retval.data = retval.data >= arg
            return retval
        elif isinstance(arg,nddata):
            retval,B = self.aligndata(arg)
            retval.data = retval.data >= B.data
            return retval
        elif isscalar(arg):
            retval = self.copy()
            retval.data = retval.data >= arg
            return retval
        else:
            raise ValueError("I don't know what to do with an argument of type"+repr(type(arg)))
    def __mul__(self,arg):
        #{{{ do scalar multiplication
        if isscalar(arg):
            #print "multiplying",self.data.dtype,"with scalar of type",type(arg)
            A = self.copy()
            if type(arg) is complex and self.data.dtype not in [complex128,complex64]:
                A.data = complex128(A.data)
            A.data *= arg
            if A.get_error() != None:
                error = A.get_error()
                error *= abs(arg)
            return A
        #}}}
        #{{{ shape and multiply
        try:
            A,B = self.aligndata(arg)
        except:
            raise CustomError("Error aligning right (arg)",arg.name(),"with left (self)",self.name())
        retval = A.copy()
        retval.data = A.data * B.data
        #}}}
        #{{{ if we have error for both the sets of data, I should propagate that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0 # we can have error on one or both, so we're going to need to add up the variances
        if Aerr != None:
            Rerr += (Aerr * B.data)**2
        if Berr != None:
            Rerr += (Berr * A.data)**2
        Rerr = sqrt(real(Rerr)) # convert back to stdev
        if Aerr == None and Berr == None:
            Rerr = None
        #}}}
        retval.set_error(Rerr)
        return retval
    def __pow__(self,arg):
        if arg == -1:
            x = self.get_error()
            result = self.copy()
            result.data = 1.0/result.data
            if x != None:
                result.set_error(abs(x.copy()/(self.data**2)))
            return result
        else:
            if self.get_error() != None:
                raise CustomError("nothing but -1 supported yet!")
            else:
                result = self.copy()
                result.data = result.data**arg
                return result
    def __div__(self,arg):
        if isscalar(arg):
            A = self.copy()
            A.data /= arg
            if A.get_error() != None:
                error = A.get_error()
                error /= abs(arg)
            return A
        try:
            A,B = self.aligndata(arg)
        except:
            raise CustomError("Error aligning right (arg) name:",arg.name(),"with left (self) name:",self.name(),"shapes are (resp):",ndshape(arg),ndshape(self))
        retval = A.copy()
        retval.data = A.data / B.data
        #{{{ if we have error for both the sets of data, I should propagate that error
        Aerr = A.get_error()
        Berr = B.get_error()
        Rerr = 0.0 # we can have error on one or both, so we're going to need to add up the variances
        dt128 = dtype('complex128')
        if Aerr != None:
            if (A.data.dtype is dt128) or (B.data.dtype is dt128):# this should avoid the error that Ryan gets
                Rerr += (complex128(Aerr)/complex128(B.data))**2
            else:
                Rerr += (Aerr/B.data)**2
        if Berr != None:
            if (A.data.dtype is dt128) or (Berr.dtype is dt128) or (B.data.dtype is dt128):# this should avoid the error that Ryan gets
                Rerr += (complex128(A.data)*complex128(Berr)/(complex128(B.data)**2))**2
            else:
                try:
                    Rerr += (A.data*Berr/(B.data**2))**2
                except:
                    raise CustomError('self was',self,
                            'arg was',arg,
                            'dtype of A.data',A.data.dtype,
                            'dtype of Berr',Berr.dtype,
                            'dtype of B.data',Berr)
        try:
            Rerr = sqrt(real(Rerr)) # convert back to stdev --> note that this has problems with complex numbers, hence the "abs" above
        except:
            raise CustomError("Rerr gave an attribute error when you passed",Rerr)
        #print "Rerr dtype",Rerr.dtype
        if Aerr == None and Berr == None:
            Rerr = None
        #}}}
        retval.set_error(Rerr)
        return retval
    def __invert__(self):
        if self.data.dtype is dtype('bool'):
            self.data = ~self.data
            return self
        else:
            raise ValueError('invert only implemented for boolean now')
    def __abs__(self):
        return self.runcopy(abs)
    __radd__ = __add__
    __rmul__ = __mul__
    def __rsub__(self,arg):
        return -1*(self-arg)
    def __neg__(self):
        return -1*self
    def __rdiv__(self,arg):
        return arg * (self**(-1))
    def real(self):
        self.data = real(self.data)
        return self
    #}}}

    #{{{ align data
    def aligndata(self,arg):
        r'''This now just returns selfout,argout
        which are aligned to each other, and which contain
        axis labels and axis errors for both'''
        #{{{ if zero dimensional, fake a singleton dimension and recurse
        #{{{ unless both are zero dimensional, in which case, just leave alone
        if ndshape(self).zero_dimensional and ndshape(arg).zero_dimensional:
            return self.copy(),arg.copy()
        #}}}
        elif ndshape(self).zero_dimensional:
            A = self.copy()
            A.dimlabels = [arg.dimlabels[0]]
            A.data = A.data.reshape(1)
            return A.aligndata(arg)
        elif ndshape(arg).zero_dimensional:
            arg = arg.copy()
            arg.dimlabels = [self.dimlabels[0]]
            arg.data = arg.data.reshape(1)
            return self.aligndata(arg)
        #}}}
        augmentdims = [x for x in arg.dimlabels if x in set(self.dimlabels)^set(arg.dimlabels)] # dims in arg but now self, ordered as they were in arg
        newdims = self.dimlabels + augmentdims # this should return the new dimensions with the order of self preserved, followed by augmentdims
        selfout = self.copy() # copy self
        if len(selfout.data.shape) == 0:
            if len(selfout.dimlabels) == 1:
                selfout.data.reshape(1)
            else:
                raise ValueError("This instance is zero dimensional (It's %s)!!!"%repr(self))
        selfshape = list(selfout.data.shape)+list(ones(len(augmentdims))) # there is no need to transpose self, since its order is preserved
        new_arg_labels = [x for x in newdims if x in arg.dimlabels] # only the labels valid for B, ordered as they are in newdims
        argout = arg.copy()
        if len(argout.data.shape) == 0:
            if len(argout.dimlabels) == 1:
                argout.data = argout.data.reshape(1)
                if argout.get_error() is not None:
                    try:
                        argout.set_error(argout.get_error().reshape(1))
                    except:
                        raise ValueError("error was"+repr(argout.get_error()))
                #print "DEBUG: I reshaped argout, so that it now has shape",argout.data.shape,"dimlabels",argout.dimlabels,"and ndshape",ndshape(argout)
            else:
                raise ValueError("The argument is zero dimensional (It's %s), but the self is not (It's %s)!!!"%(repr(arg),repr(self)))
        argshape = list(ones(len(newdims)))
        #print "DEBUG 2: shape of self",ndshape(self),"self data shape",self.data.shape,"shape of arg",ndshape(arg),"arg data shape",arg.data.shape
        #print "DEBUG 3: shape of selfout",ndshape(selfout),"selfout data shape",selfout.data.shape,"shape of argout",ndshape(argout),"argout data shape",argout.data.shape
        #{{{ wherever the dimension exists in arg, pull the shape from arg
        for j,k in enumerate(newdims):
            if k in argout.dimlabels:
                try:
                    argshape[j] = argout.data.shape[argout.axn(k)]
                except:
                    raise CustomError("There seems to be a problem because the shape of argout is now len:%d "%len(argout.data.shape),argout.data.shape,"while the dimlabels is len:%d "%len(argout.dimlabels),argout.dimlabels)
        #}}}
        argorder = map(argout.dimlabels.index,new_arg_labels) # for each new dimension, determine the position of the original dimension
        selfout.data = selfout.data.reshape(selfshape) # and reshape to its new shape
        selfout.dimlabels = newdims
        try:
            argout.data = argout.data.transpose(argorder).reshape(argshape) # and reshape the data
        except ValueError,Argument:
            raise ValueError('the shape of the data is '+repr(argout.data.shape)+' the transpose '+repr(argorder)+' and the new shape '+repr(argshape)+' original arg: '+repr(Argument))
        argout.dimlabels = newdims
        if selfout.get_error() != None:
            try:
                temp = selfout.get_error().copy().reshape(selfshape)
            except ValueError,Argument:
                raise ValueError("The instance (selfout) has a shape of "+repr(selfout.data.shape)+" but its error has a shape of"+repr(selfout.get_error().shape)+"!!!\n\n(original argument:\n"+repr(Argument)+"\n)")
            selfout.set_error(temp)
        if argout.get_error() != None:
            try:
                temp = argout.get_error().copy().transpose(argorder).reshape(argshape)
            except ValueError,Argument:
                raise ValueError("The argument (argout) has a shape of "+repr(argout.data.shape)+" but its error has a shape of"+repr(argout.get_error().shape)+"(it's "+repr(argout.get_error())+")!!!\n\n(original argument:\n"+repr(Argument)+"\n)")
            argout.set_error(temp)
        if (len(selfout.axis_coords)>0) or (len(argout.axis_coords)>0):
            #{{{ transfer the errors and the axis labels
            #{{{ make dictionaries for both, and update with info from both, giving preference to self
            axesdict = selfout.mkd()
            errordict = selfout.mkd()
            #{{{ add the axes and errors for B
            if type(arg.axis_coords) is list:
                if len(arg.axis_coords) > 0:
                    axesdict.update(arg.mkd(arg.axis_coords))
            if type(arg.axis_coords_error) is list:
                if len(arg.axis_coords_error) > 0 and not all([x is None for x in arg.axis_coords_error]):
                    errordict.update(arg.mkd(arg.axis_coords_error))
            #}}}
            #{{{ add the errors for A
            if type(self.axis_coords) is list:
                if len(self.axis_coords) > 0:
                    axesdict.update(self.mkd(self.axis_coords))
            if type(self.axis_coords_error) is list:
                if len(self.axis_coords_error) > 0 and not all([x is None for x in self.axis_coords_error]):
                    errordict.update(self.mkd(self.axis_coords_error))
            #}}}
            #}}}
            selfout.axis_coords_error = selfout.fld(errordict)
            argout.axis_coords_error = selfout.fld(errordict)
            selfout.axis_coords = selfout.fld(axesdict)
            argout.axis_coords = selfout.fld(axesdict)
            #}}}
            selfout.axis_coords_units = [None]*len(newdims)
            argout.axis_coords_units = [None]*len(newdims)
            for thisdim in newdims:
                if thisdim in self.dimlabels:
                    selfout.set_units(thisdim,self.get_units(thisdim))
                    argout.set_units(thisdim,self.get_units(thisdim))
                elif thisdim in arg.dimlabels:
                    selfout.set_units(thisdim,arg.get_units(thisdim))
                    argout.set_units(thisdim,arg.get_units(thisdim))
        return selfout,argout
    #}}}

    def removeDuplicates(self,thisAxis):
        """ Remove the duplicates along a given dimension. Averages the data values corresponding to the duplicates, sets error to the std of the duplicate entries.
        Limitations:
        Data must be numeric, cannot be str
        Only handles 1-d data sets.

        thisaxis - str - axis to remove duplicates along 

        returns new nddata must call as newSet = oldSet.removeDuplicates('axisName')

        """
        indepAxis = list(set(self.getaxis(thisAxis)))
        dataDim = []
        errorDim = []
        for count,axisVal in enumerate(indepAxis):
            dataVals = self[thisAxis,lambda x: x==axisVal].data
            if len(dataVals) > 1:
                dataDim.append(average(dataVals))
                errorDim.append(std(dataVals))
            else:
                dataDim.append(float(dataVals))
                if self.get_error() is not None:
                    errorDim.append(float(self[thisAxis,lambda x: x==axisVal].get_error()))
        if len(errorDim) > 1:
            selfout = nddata(array(dataDim)).rename('value',thisAxis).labels(thisAxis,array(indepAxis)).set_error(array(errorDim)) 
        else:
            selfout = nddata(array(dataDim)).rename('value',thisAxis).labels(thisAxis,array(indepAxis))
        return selfout


    #{{{ integrate, differentiate, and sum
    def integrate(self,thisaxis,backwards = False):
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        self.run_nopop(cumsum,thisaxis)
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        if len(self.axis_coords)>0:
            t = self.getaxis(thisaxis)
            dt = t[1]-t[0]
            self.data *= dt
        return self

    def diff(self,thisaxis,backwards = False):
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        self.run_nopop(mydiff,thisaxis)
        if backwards is True:
            self.data = self[thisaxis,::-1].data
        if len(self.axis_coords)>0:
            t = self.getaxis(thisaxis)
            dt = t[1]-t[0]
            self.data /= dt
        return self

    def sum(self,axes):
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print '|-ERROR FINDING DIMENSION-----'
                print '| dimlabels is: ',self.dimlabels
                print "| doesn't contain: ",axes[j]
                print '|-----------------------------'
                raise
            self.data = sum(self.data,
                    axis=thisindex)
            self._pop_axis_info(thisindex)
        return self

    def sum_nopop(self,axes):
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            temp = list(self.data.shape)
            temp[thisindex] = 1
            self.data = sum(self.data,
                    axis=thisindex)
            self.data = self.data.reshape(temp)
        return self
    #}}}


    #{{{ max and mean
    def _wrapaxisfuncs(self,func):
        #{{{ for convenience, wrap the max and min functions
        if func == max:
            func = amax
        if func == min:
            func = amin
        if func == diff:
            func = mydiff
        return func
        #}}}

    def argmax(self,*args,**kwargs):
        r'find the max along a particular axis, and get rid of that axis, replacing it with the index number of the max value'
        #{{{ process arguments
        axes = self._possibly_one_axis(*args)
        raw_index = False
        if 'raw_index' in kwargs.keys():
            raw_index = kwargs.pop('raw_index')
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        if (type(axes) is str):
            axes = [axes]
        #}}}
        for j in range(0,len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            if raw_index:
                self.data = argmax(self.data,
                        axis=thisindex)
            else:
                self.data = self.axis_coords[thisindex][argmax(self.data,
                    axis=thisindex)]
            self._pop_axis_info(thisindex)
        return self

    def argmin(self,axes,raw_index = False):
        r'find the min along a particular axis, and get rid of that axis, replacing it with the index number of the min value'
        if (type(axes) is str):
            axes = [axes]
        for j in range(0,len(axes)):
            try:
                thisindex = self.axn(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            if raw_index:
                self.data = argmin(self.data,
                        axis=thisindex)
            else:
                self.data = self.axis_coords[thisindex][argmin(self.data,
                    axis=thisindex)]
            self._pop_axis_info(thisindex)
        return self

    def mean_all_but(self,listofdims):
        'take the mean over all dimensions not in the list'
        for dimname in self.dimlabels:
            if not (dimname in listofdims):
                self.mean(dimname)
        return self

    def mean(self,*args,**kwargs):
        r'Take the mean and set the error to the standard deviation'
        #{{{ process arguments
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        axes = self._possibly_one_axis(*args)
        return_error = False
        if return_error in kwargs.keys():
            return_error = kwargs.pop('return_error')
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        if (type(axes) is str):
            axes = [axes]
        #}}}
        for j in range(0,len(axes)):
            try:
                thisindex = self.dimlabels.index(axes[j])
            except:
                print 'error, dimlabels is: ',self.dimlabels
                print "doesn't contain: ",axes[j]
                raise
            if self.data_error is not None:
                this_axis_length = self.data.shape[thisindex]
                try:
                    self.data_error = sqrt(sum((self.data*self.data_error)**2,
                            axis=thisindex)/(this_axis_length**2))
                except:
                    raise CustomError('shape of data',shape(self.data),'shape of data error',shape(self.data_error))
            if return_error: # since I think this is causing an error
                thiserror = std(self.data,
                        axis=thisindex)
                if isscalar(thiserror):
                    thiserror = r_[thiserror]
                self.set_error(thiserror) # set the error to the standard deviation
            self.data = mean(self.data,
                    axis=thisindex)
            self._pop_axis_info(thisindex)
        return self

    def mean_nopop(self,axis):
        self = self.run_nopop(mean,axis=axis)
        return self
    #}}}

    #{{{ running functions and popping dimensions
    def _pop_axis_info(self,thisindex):
        r'pop axis by index'
        self.dimlabels.pop(thisindex)
        if self.axis_coords!=[]:
            self.axis_coords.pop(thisindex)
            try:
                self.axis_coords_error.pop(thisindex)
            except:
                raise CustomError('trying to pop',thisindex,'from',self.axis_coords_error)
            if len(self.axis_coords_units) > 0:
                try:
                    self.axis_coords_units.pop(thisindex)
                except:
                    raise CustomError('trying to pop',thisindex,'from',self.axis_coords_units)
        return self

    def popdim(self,dimname):
        thisindex = self.axn(dimname)
        thisshape = list(self.data.shape)
        if thisshape[thisindex]!=1:
            raise CustomError("trying to pop a dim that's not length 1")
        thisshape.pop(thisindex)
        self.data = self.data.reshape(thisshape)
        self._pop_axis_info(thisindex)
        return self

    def cropped_log(self,magnitude = 4):
        r'''For the purposes of plotting, this generates a copy where I take the log, spanning "magnitude" orders of magnitude
        This is designed to be called as abs(instance).cropped_log(), so it doesn't make a copy'''
        self = self.run(log10)
        self -= self.data.flatten().max() - 4 # span only 4 orders of magnitude
        self.data[self.data < 0] = 0
        return self

    def runcopy(self,*args):
        newdata = self.copy()
        func = args[0]
        func = self._wrapaxisfuncs(func)
        if len(args)>1:
            axis = args[1]
            thisindex = newdata.dimlabels.index(axis)
            newdata.data = func(newdata.data,axis=thisindex)
            newdata._pop_axis_info(thisindex)
        else:
            newdata.data = func(newdata.data)
        return newdata

    def run(self,*args):
        func = args[0]
        func = self._wrapaxisfuncs(func)
        if len(args)>1:
            axis = args[1]
            try:
                thisindex = self.dimlabels.index(axis)
            except:
                if type(axis) is not str:
                    raise ValueError('The format of run is run(func,"axisname"), but you didn\'t give a string as the second argument -- maybe you fed the arguments backwards?')
            self.data = func(self.data,axis=thisindex)
            self._pop_axis_info(thisindex)
        else:
            self.data = func(self.data)
        return self

    def run_nopop(self,func,axis):
        func = self._wrapaxisfuncs(func)
        try:
            thisaxis = self.dimlabels.index(axis)
        except:
            raise CustomError("I couldn't find the dimension",axis,"in the list of axes",self.dimlabels)
        temp = list(self.data.shape)
        temp[thisaxis] = 1
        numnonoptargs = len(getargspec(func)[0])-len(getargspec(func)[3])
        if numnonoptargs == 1:
            try:
                self.data = func(self.data,axis=thisaxis)
            except TypeError:
                self.data = func(self.data,axes=thisaxis)
        elif numnonoptargs == 2:
            try:
                self.data = func(self.getaxis(axis),self.data,axis=thisaxis)
            except TypeError:
                self.data = func(self.getaxis(axis),self.data,axes=thisaxis)
        else:
            raise CustomError('you passed a function to run_nopop that doesn\'t have either one or two arguments!')
        #{{{ if the function doesn't rip out the dim, make sure we don't change the dims
        if len(self.data.shape)==len(temp):
            temp[thisaxis] = self.data.shape[thisaxis]
        #}}}
        self.data = self.data.reshape(temp)
        return self
    #}}}

    #{{{ ft-related functions
    def convolve(self,axisname,filterwidth,convfunc = (lambda x,y: exp(-(x**2)/(2.0*(y**2))))):
        r'''perform a normalized convolution'''
        #{{{ make a version of x that is oriented along the correct dimension
        x = self.getaxis(axisname).copy()
        x_centerpoint = (x[-1]+x[0])/2
        x -= x_centerpoint # so that zero is in the center
        x = ifftshift(x) # so that it's convolved about time 0
        thisaxis = self.axn(axisname)
        #}}}
        myfilter = convfunc(x,filterwidth)
        myfilter /= myfilter.sum()
        filtershape = ones_like(self.data.shape)
        filtershape[thisaxis] = len(myfilter)
        myfilter = myfilter.reshape(filtershape)
        self.data = ifft(fft(self.data,axis = thisaxis)*fft(myfilter,axis = thisaxis),axis = thisaxis)
        return self

    def _ft_conj(self,x):
        pairs = [('s','Hz'),('m',r'm^{-1}')]
        a,b = zip(*tuple(pairs))
        if x in a:
            return b[a.index(x)]
        elif x in b:
            return a[b.index(x)]
        else:
            return None

    def ftshift(self,axis):
        self.data = fftshift(self.data,axes = self.axn(axis))
        x = self.getaxis(axis)
        x[:] = fftshift(x)
        j = len(x)/2 # given the floor, this works out to be the central index
        x_subset = x[:j]
        x_subset -= x_subset[-1] + x[j+1] # zero and set to this
        return self

    def ft(self,*args,**kwargs):
        self.set_prop('FT',True)
        #{{{ process arguments
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        axes = self._possibly_one_axis(*args)
        #kwargs: shiftornot=False,shift=None,pad = False
        shiftornot,shift,pad,automix = self._process_kwargs([('shiftornot',False),('shift',None),('pad',False),('automix',False)],**kwargs)
        if shift != None:
            shiftornot = shift
        if (type(axes) is str):
            axes = [axes]
        if not (type(shiftornot) is list):
            shiftornot = [bool(shiftornot)]*len(axes)
        #}}}
        for j in range(0,len(axes)):
            if self.get_units(axes[j]) is not None:
                self.set_units(axes[j],self._ft_conj(self.get_units(axes[j])))
            try:
                thisaxis = self.dimlabels.index(axes[j])
            except:
                raise CustomError('error, dimlabels is: ',self.dimlabels)
            padded_length = self.data.shape[thisaxis]
            if pad is True:
                padded_length = 2**(ceil(log2(padded_length)))
            elif pad:
                padded_length = pad
            self.data = fft(self.data,n = padded_length,axis=thisaxis)
            if bool(shiftornot[j]):
                if automix:
                    raise ValueError("You can't use automix and shift at the same time --> it doesn't make sense")
                self.data = fftshift(self.data,axes=[thisaxis])
            t = self.getaxis(axes[j])
            if t is not None:
                dt = t[1]-t[0] # the dwell gives the bandwidth, whether or not it has been zero padded
                self.ft_start_time = t[0]
                self.data *= dt
                self.axis_coords[thisaxis] = linspace(0,1./dt,padded_length)
                if bool(shiftornot[j]):
                    mask = self.axis_coords[thisaxis] > 0.5/dt
                    #{{{ just the axis part of ftshift
                    x = self.axis_coords[thisaxis]
                    x[:] = fftshift(x)
                    j = len(x)/2 # given the floor, this works out to be the central index
                    x_subset = x[:j]
                    x_subset -= x_subset[-1] + x[j+1] # zero and set to this
                    #}}}
            if automix:
                sw = 1.0/dt
                carrier = abs(self).mean_all_but(axes[j]).argmax(axes[j]).data
                print "I find carrier at",carrier
                add_to_axis = (automix - carrier) / sw
                print "which is",add_to_axis,"times the sw of",sw,"off from the automix value of",automix
                x = self.getaxis(axes[j])
                x += round(add_to_axis)*sw
        return self

    def ift(self,*args,**kwargs):
        self.set_prop('FT',False)
        #{{{ process arguments
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        axes = self._possibly_one_axis(*args)
        #kwargs: shiftornot=False,shift=None,pad = False
        shiftornot,shift,pad = self._process_kwargs([('shiftornot',False),('shift',None),('pad',False)],**kwargs)
        if shift != None:
            shiftornot = shift
        if (type(axes) is str):
            axes = [axes]
        if not (type(shiftornot) is list):
            shiftornot = [bool(shiftornot)]*len(axes)
        #}}}
        for j in range(0,len(axes)):
            if self.get_units(axes[j]) is not None:
                self.set_units(axes[j],self._ft_conj(self.get_units(axes[j])))
            try:
                thisaxis = self.dimlabels.index(axes[j])
            except:
                raise CustomError('error, dimlabels is: ',self.dimlabels)
            padded_length = self.data.shape[thisaxis]
            if pad is True:
                padded_length = 2**(ceil(log2(padded_length)))
            elif pad:
                padded_length = pad
            if bool(shiftornot[j]):
                self.data = ifftshift(self.data,axes=[thisaxis])
            self.data = ifft(self.data,n = padded_length,axis=thisaxis)
            t = self.getaxis(axes[j])
            if t is not None:
                dt = t[1]-t[0]
                self.data *= size(t) * dt # here, the algorithm divides by N, so for integration, we need to not do that
                #{{{ shiftornot specifies the shifting of the initial ft, not this result, so we always return a 0->1 time axis
                self.axis_coords[thisaxis] = linspace(0,1./dt,size(t)) + self.ft_start_time # note that I offset by ft_start_time, which I pull from when I ft'd
                #}}}
        return self
    #}}}

    #{{{ interpolation and binning
    def run_avg(self,thisaxisname,decimation = 20,centered = False):
        'a simple running average'
        temp = self.getaxis(thisaxisname).size % decimation
        decimation = int(decimation)
        if temp != 0:
            if centered:
                self = self[thisaxisname,temp/2:-int(temp/2. + 0.5)]
            else:
                self = self[thisaxisname,0:-temp]
        thisaxis = nddata(self.getaxis(thisaxisname),[-1],[thisaxisname])
        self.setaxis(thisaxisname,[])
        self.chunkoff(thisaxisname,['avg'],[decimation])
        self.run(mean,'avg')
        thisaxis.chunkoff(thisaxisname,['avg'],[decimation])
        thisaxis.run(mean,'avg')
        self.setaxis(thisaxisname,thisaxis.data)
        return self

    def histogram(*args,**kwargs):
        per_hit = 1e-3
        if 'per_hit' in kwargs.keys():
            per_hit = kwargs.pop('per_hit')
        if len(kwargs) > 0:
            raise ValueError("I don't understand the arguments:"+repr(kwargs))
        if len(args) == 1:
            if type(args[0]) is ndarray:
                if args[0].shape[2] != len(self.dimlabels):
                    raise ValueError("You must pass an N x M array, where M is the number of dimensions in this array!")
            elif type(args[0]) is dict:
                raise ValueError("should eventually support dictionaries (via fld), but doesn't yet")
            else:
                raise ValueError("I don't know what funny business you're up to passing me a"+repr(type(args[0])))
        else:
            raise ValueError("should eventually support array, label pair, but doesn't yet")

    def interp(self,axis,axisvalues,past_bounds = None,verbose = False,**kwargs):
        'interpolate data values given axis values'
        oldaxis = self.getaxis(axis)
        if (type(axisvalues) is int) or (type(axisvalues) is int32):
            axisvalues = linspace(oldaxis[0],oldaxis[-1],axisvalues)
        elif isscalar(axisvalues):
            axisvalues = r_[axisvalues]
        elif (type(axisvalues) not in [ndarray,tuple]):
            raise ValueError("You passed a target axis of type"+repr(type(axisvalues))+"which I don't get")
        if any(imag(axisvalues) > 1e-38):
            raise ValueError("I can't interpolate imaginary values")
        else:
            axisvalues = real(axisvalues)
        if past_bounds == None:
            axisvalues[axisvalues<oldaxis.min()] = oldaxis.min()
            axisvalues[axisvalues>oldaxis.max()] = oldaxis.max()
        elif not (past_bounds == 'fail'):
            if type(past_bounds) is tuple:
                if len(past_bounds) == 2:
                    axisvalues[axisvalues<oldaxis.min()] = past_bounds[0]
                    axisvalues[axisvalues>oldaxis.max()] = past_bounds[1]
                else:
                    raise TypeError('If you pass axisvalues as a tuple, it must be of length 2!')
            else:
                axisvalues[axisvalues<oldaxis.min()] = past_bounds
                axisvalues[axisvalues>oldaxis.max()] = past_bounds
        rdata = real(self.data)
        idata = imag(self.data)
        if 'kind' in kwargs.keys():
            thiskind = kwargs.pop('kind')
        else:
            thiskind = 'cubic'
            if len(rdata) < 4:
                thiskind = 'quadratic'
                if len(rdata) < 3:
                    thiskind = 'linear'
        thisaxis = self.axn(axis)
        if verbose: print 'Using %s interpolation'%thiskind
        interpfunc =  interp1d(oldaxis,rdata,kind = thiskind,axis = thisaxis)
        try:
            rdata = interpfunc(axisvalues)
        except:
            raise CustomError("dtype of axis is"+repr(axisvalues.dtype))
        interpfunc =  interp1d(oldaxis,idata,kind = thiskind,axis = thisaxis)
        idata = interpfunc(axisvalues)
        self.data = rdata + 1j * idata
        self.setaxis(axis,axisvalues)
        return self

    def invinterp(self,axis,values,**kwargs):
        'interpolate axis values given data values'
        copy = False
        if 'copy' in kwargs.keys():
            copy = kwargs.pop('copy')
        if isscalar(values):
            values = r_[values]
        origdata = self.data.copy()
        origaxis = self.getaxis(axis).copy()
        if any(imag(values) > 1e-38):
            raise ValueError("I can't interpolate imaginary values")
        else:
            values = real(values)
        args = origdata.argsort()
        origdata = origdata[args]
        rdata = real(origaxis)
        idata = imag(origaxis)
        rdata = rdata[args]
        idata = idata[args]
        #{{{ determine type of interpolation
        if 'kind' in kwargs.keys():
            thiskind = kwargs.pop('kind')
        else:
            thiskind = 'cubic'
            if len(rdata) < 4:
                thiskind = 'quadratic'
                if len(rdata) < 3:
                    thiskind = 'linear'
        #}}}
        interpfunc =  interp1d(origdata,rdata,kind = thiskind)
        try:
            rdata = interpfunc(values)
        except:
            raise CustomError('You passed',values,'and the data spans from',origdata.min(),'to',origdata.max())
        interpfunc =  interp1d(origdata,idata,kind = thiskind)
        idata = interpfunc(values)
        cdata = rdata + 1j * idata
        if copy:
            newaxis = 'name data before interpolation'
            if self.name() is not None:
                newaxis = self.name()
            retval = nddata(cdata,[-1],[newaxis])
            retval.labels(newaxis,values)
            return retval
        else:
            self.data = values
            self.setaxis(axis,cdata)
            return self
    #}}}

    def multimin(self,minfunc,axisname,filterwidth,numberofmins):
        cost = self.copy().convolve(axisname,filterwidth).run_nopop(minfunc)
        for j in range(0,numberofmins):
            #{{{ find the x value at which the minimum occurs
            xvalues = cost.copy().argmin(axisname,raw_index = True)
            #}}}

    def repwlabels(self,axis):
        return None

    #{{{ functions to manipulate and return the axes
    def reorder(self,*axes,**kwargs):
        first = True
        if 'first' in kwargs:
            first = kwargs.pop('first')
        if len(kwargs) > 0:
            raise ValueError("I don't understand your kwargs!")
        if len(axes) == 1:
            axes = axes[0]
        else:
            axes = axes
        if type(axes) is str:
            axes = [axes]
        if len(axes) < len(self.dimlabels):
            if first:
                oldorder = list(self.dimlabels)
                for thisaxis in axes:
                    oldorder.pop(oldorder.index(thisaxis))
                axes = axes + oldorder
            else:
                raise ValueError("False (to put these axes at the end) is not yet supported")
        try:
            neworder = map(self.dimlabels.index,axes)
        except ValueError:
            raise CustomError('one of',axes,'not in',self.dimlabels)
        self.dimlabels = map(self.dimlabels.__getitem__,neworder)
        if len(self.axis_coords)>0:
            try:
                self.axis_coords = map(self.axis_coords.__getitem__,neworder)
            except:
                raise CustomError('problem mapping',map(len,self.axis_coords),'onto',neworder)
            if len(self.axis_coords_units)>0:
                try:
                    self.axis_coords_units = map(self.axis_coords_units.__getitem__,neworder)
                except:
                    raise CustomError('problem mapping',map(len,self.axis_coords_units),'onto',neworder)
        try:
            self.data = self.data.transpose(neworder)
        except ValueError:
            raise CustomError('you can\'t reorder',self.dimlabels,'as',neworder)
        if self.data_error is not None:
            self.data_error = self.data_error.transpose(neworder)
        return self
    def plot_labels(self,labels,fmt = None,**kwargs_passed):
        r'this only works for one axis now'
        axisname = self.dimlabels[0]
        if fmt is None:
            plot_label_points(self.getaxis(axisname),self.data,labels,**kwargs_passed)
        else:
            plot_label_points(self.getaxis(axisname),self.data,[fmt%j for j in labels],**kwargs_passed)
        return
    def labels(self,*args):
        r'''label the dimensions, given in listofstrings with the axis labels given in listofaxes -- listofaxes must be a numpy array;
        you can pass either a list or a axis name (string)/axis label (numpy array) pair'''
        if len(args) == 2:
            listofstrings,listofaxes = args
        elif len(args) == 1 and type(args[0]) is dict:
            listofstrings = args[0].keys()
            listofaxes = args[0].values()
        else:
            raise CustomError("I can't figure out how to deal with the arguments",args)
        for j in range(0,len(listofaxes)):
            if type(listofaxes[j]) is list:
                listofaxes[j] = array(listofaxes[j])
        if len(self.axis_coords) == 0:
            self.axis_coords = [[]]*len(self.dimlabels)
            self.axis_coords_error = [None]*len(self.dimlabels)
        if type(listofstrings) is str:
            listofstrings = [listofstrings]
            listofaxes = [listofaxes]
        if type(listofstrings) is not list:
            raise TypeError("the arguments passed to the .labels() method must be a list of the axis names followed by the list of the axis arrays")
        elif all(map(( lambda x: type(x) is str_ ),listofstrings)):
            listofstrings = map(str,listofstrings)
        elif not all(map(( lambda x: type(x) in [str,str_] ),listofstrings)):
            raise TypeError("the arguments passed to the .labels() method must be a list of the axis names followed by the list of the axis arrays")
        for j in range(0,len(listofstrings)):
            if listofaxes[j] is None:
                self.axis_coords[self.axn(listofstrings[j])] = None
            else:
                #{{{ test that the axis is the right size
                if isscalar(listofaxes[j]):#interpret as a timestep
                    listofaxes[j] = listofaxes[j] * r_[0:ndshape(self)[listofstrings[j]]]
                if type(listofaxes[j]) not in [ndarray,list]:
                    raise TypeError('You passed an axis label of type '+repr(type(listofaxes[j]))+' for the axis '+listofstrings[j]+' to the labels method, which you can\'t do --> it must be an nddata')
                if (len(listofaxes[j]) != ndshape(self)[listofstrings[j]]) and (len(listofaxes[j])!=0):
                    raise IndexError("You're trying to attach an axis of len %d to the '%s' dimension, which has %d data points"%(len(listofaxes[j]),listofstrings[j],ndshape(self)[listofstrings[j]]))
                #}}}
                try:
                    self.axis_coords[self.dimlabels.index(listofstrings[j])] = listofaxes[j]
                except:
                    try:
                        raise CustomError("Can't assign the coordinates to "+listofstrings[j]+" as "+listofaxes[j].__repr__())
                    except:
                        raise CustomError("length of listofaxes (",len(listofaxes),") isn't same length as ",listofstrings)
        return self
    def check_axis_coords_errors(self):
        if len(self.axis_coords_error) > len(self.dimlabels):
            raise ValueError("this failed because there are more sets of axis errors than there are axes!\nlen(axis_coords_error) = %s\naxes = %s"%(repr(len(self.axis_coords_error)),repr(self.dimlabels)))
    def sort(self,axisname):
        whichaxis = self.dimlabels.index(axisname)
        order = argsort(self.axis_coords[whichaxis])
        datacopy = self.copy()
        for j in range(0,len(order)): # do it this way, so that it deals with other dimensions correctly
            self.check_axis_coords_errors()
            self[axisname,j] = datacopy[axisname,order[j]]
        self.axis_coords[whichaxis] = self.axis_coords[whichaxis][order]
        return self
    def copyaxes(self,other):
        # in the case that the dimensions match, and we want to copy the labels
        self.axis_coords = other.axis_coords
        return self
    def axis(self,axisname):
        'returns a 1-D axis for further manipulation'
        thisaxis = self.axn(axisname)
        return nddata(self.getaxis(axisname).copy(),[-1],[axisname]).labels(axisname,self.getaxis(axisname).copy())
    def _axis_inshape(self,axisname):
        newshape = ones(len(self.data.shape),dtype = 'uint')
        thisaxis = self.axn(axisname)
        newshape[thisaxis] = self.data.shape[thisaxis]
        newshape = list(newshape)
        return self.getaxis(axisname).copy().reshape(newshape)
    def retaxis(self,axisname):
        thisaxis = self._axis_inshape(str(axisname))
        return nddata(thisaxis,thisaxis.shape,list(self.dimlabels)).labels(str(axisname),thisaxis.flatten())
    def fromaxis(self,*args,**kwargs):
        '''enter just the axis, to return the axis,
        or enter a list of axisnames, followed by a function to act on them'''
        overwrite = False
        if 'overwrite' in kwargs.keys() and kwargs['overwrite'] == True:
            overwrite = True
        as_array = False
        if 'as_array' in kwargs.keys() and kwargs['as_array'] == True:
            as_array = True
        if len(args) == 1:
            if type(args[0]) is str:
                retval = self.retaxis(args[0])
                #{{{ copied from old retaxis function, then added the overwrite capability
                thisaxis = self._axis_inshape(axisname)
                if overwrite:
                    self.data = thisaxis
                else:
                    retval = nddata(thisaxis,thisaxis.shape,list(self.dimlabels)).labels(axisname,thisaxis.flatten())
                    retval.axis_coords_units = list(self.axis_coords_units)
                    retval.data_units = self.data_units
                    return retval
                #}}}
            else:
                raise CustomError("I don't know what to do!")
        if len(args) == 2:
            axisnames = args[0]
            func = args[1]
            if type(axisnames) is not list:
                axisnames = [axisnames]
        else:
            raise CustomError('Wrong number of arguments!!')
        if issympy(func):
            func = sympy.lambdify(*tuple(map(sympy.var,axisnames) + [func,"numpy"]))
        if func.func_code.co_argcount != len(axisnames):
            raise CustomError("The axisnames you passed and the argument count don't match")
        list_of_axes = [self._axis_inshape(x) for x in axisnames]
        retval = func(*list_of_axes)
        newshape = ones_like(list_of_axes[0].shape)
        for j in list_of_axes:
            newshape *= array(j.shape)
        if overwrite:
            self.data = retval.reshape(newshape)
            return self
        else:
            if as_array:
                return retval.reshape(newshape)
            else:
                retval =  nddata(retval,
                        newshape,
                        list(self.dimlabels)).labels(axisnames,
                                [self.getaxis(x).copy() for x in axisnames])
                retval.axis_coords_units = list(self.axis_coords_units)
                retval.data_units = self.data_units
                return retval
    def getaxis(self,axisname):
        if self.axis_coords is None or len(self.axis_coords) == 0:
            return None
        else:
            retval = self.axis_coords[self.axn(axisname)]
        if retval is None:
            return None
        elif len(retval) > 0:
            return retval
        else:
            return None
    def setaxis(self,axis,value):
        if type(value) in [float,int,double,float64]:
           value = linspace(0.,value,self.axlen(axis))
        if type(value)is list:
            value = array(value)
        if ( self.axis_coords is None ) or self.axis_coords == []:
            self.axis_coords = [None] * len(self.dimlabels)
        try:
            self.axis_coords[self.axn(axis)] = value
        except:
            raise ValueError("I can't set this -- axis coords is",self.axis_coords)
    def getaxisshape(self,axisname):
        thishape = ones(len(self.dimlabels))
        thisaxis = self.dimlabels.index(axisname) 
        thishape[thisaxis] = self.data.shape[thisaxis]
        return thishape
    def circshift(self,axis,amount):
        if amount!=0:
            if abs(amount) > ndshape(self)[axis]:
                CustomError("Trying to circshift by ",amount,"which is bitter than the size of",axis)
            newdata = ndshape(self).alloc(dtype=self.data.dtype)
            newdata[axis,:-amount] = self[axis,amount:]
            newdata[axis,-amount:] = self[axis,:amount]
            self.data = newdata.data
        return self
    #}}}
    #{{{ breaking up and combining axes
    def smoosh(self,dimstocollapse,dimname = 0,verbose = False):
        r'''collapse multiple dimensions into one dimension
        if dimname is:
            None, create a new (direct product) name,
            a number, lump the existing dimension into the number given, in the list
            a string, same as the previous, where the string can be part of the list or not
        '''
        #{{{ first, put them all at the end, in order given here
        retained_dims = list(self.dimlabels)
        if verbose: print "old order",retained_dims
        #{{{ if I'm using a dimension here, be sure to grab its current position
        if dimname is None:
            final_position = -1
            dimname = ' $\\times$ '.join(dimstocollapse)
        else:
            if type(dimname) is int:
                dimname = dimstocollapse[dimname]
                final_position = self.axn(dimname)
            elif dimname in self.dimlabels:
                final_position = self.axn(dimname)
            else:
                final_position = -1
        #}}}
        for this_name in dimstocollapse:
            retained_dims.pop(retained_dims.index(this_name))
        new_order = retained_dims + dimstocollapse
        self.reorder(new_order)
        if verbose: print "new order",new_order
        #}}}
        #{{{ then, reshape the data (and error)
        if verbose: print "old shape",self.data.shape
        new_shape = list(self.data.shape)[:-len(dimstocollapse)]
        if verbose: print "dimensions to keep",new_shape
        dimstocollapse_shapes = array(self.data.shape[-len(dimstocollapse):])
        new_shape += [dimstocollapse_shapes.prod()]
        self.data = self.data.reshape(new_shape)
        if self.get_error() is not None:
            self.set_error(self.get_error().reshape(new_shape))
        if verbose: print "new shape",self.data.shape
        #}}}
        #{{{ now for the tricky part -- deal with the axis labels
        #{{{ in order, make a list of the relevant axis names, dtypes, and sizes
        axes_with_labels = [j for j in dimstocollapse if self.getaxis(j) is not None] # specifically, I am only concerned with the ones I am collapsing that have labels
        axes_with_labels_haserror = [self.get_error(j) is not None for j in axes_with_labels]
        axes_with_labels_dtype = [(j,self.getaxis(j).dtype) for j in axes_with_labels]
        axes_with_labels_size = [self.getaxis(j).size for j in axes_with_labels]
        #}}}
        if verbose: print "the dtype that I want is:",axes_with_labels_dtype
        if verbose: print "the axes that have labels are:",axes_with_labels
        if verbose: print "the axes that have labels have sizes:",axes_with_labels_size
        multidim_axis_error = None
        if len(axes_with_labels_dtype) > 0:
            # create a new axis of the appropriate shape and size
            multidim_axis_label = empty(axes_with_labels_size,dtype = axes_with_labels_dtype)
            if any(axes_with_labels_haserror):
                multidim_axis_error = empty(axes_with_labels_size,
                        dtype = [(axes_with_labels[j],self.getaxis(axes_with_labels[j]).dtype)
                            for j in range(len(axes_with_labels)) if axes_with_labels_haserror[j]])
            # one at a time index the relevant dimension, and load in the information
            full_slice = [slice(None,None,None)]*len(axes_with_labels_dtype)
            for this_index,thisdim in enumerate(axes_with_labels):
                axis_for_thisdim = self.getaxis(thisdim)
                if axes_with_labels_haserror[this_index]:
                    axis_error_for_thisdim = self.get_error(thisdim)
                if verbose: print "the axis for",thisdim,"is",axis_for_thisdim
                for j in range(axes_with_labels_size[this_index]):
                    this_slice = list(full_slice)
                    this_slice[this_index] = j # set this element
                    multidim_axis_label[thisdim][tuple(this_slice)] = axis_for_thisdim[j]
                    if axes_with_labels_haserror[this_index]:
                        multidim_axis_error[thisdim][tuple(this_slice)] = axis_error_for_thisdim[j]
            if verbose: print "shape of multidim_axis_label is now",multidim_axis_label.shape,"(",axes_with_labels,")"
            if verbose: print "multidim_axis_label is:\n",repr(multidim_axis_label)
            multidim_axis_label = multidim_axis_label.flatten() # then flatten the axis
            if verbose: print "shape of multidim_axis_label is now",multidim_axis_label.shape
            if verbose: print "multidim_axis_label is:\n",repr(multidim_axis_label)
            #{{{ create a new axis dictionary with the new info
            axis_coords_dict = self.mkd(self.axis_coords)
            axis_coords_error_dict = self.mkd(self.axis_coords_error)
            axis_coords_dict[dimname] = multidim_axis_label
            axis_coords_error_dict[dimname] = multidim_axis_error
            if verbose: print "end up with axis_coords_dict",axis_coords_dict
            if verbose: print "end up with axis_coords_error_dict",axis_coords_error_dict
            #}}}
        #}}}
        #{{{ make new dimlabels, and where relevant, project the new dictionary onto these dimlabels
        self.dimlabels = retained_dims + [dimname]
        if verbose: print "new dimlabels",self.dimlabels
        if len(axes_with_labels_dtype) > 0:
            self.axis_coords = self.fld(axis_coords_dict)
            self.axis_coords_error = self.fld(axis_coords_error_dict)
        if verbose: print "new axis coords",self.axis_coords
        if verbose: print "new axis coords errors",self.axis_coords_error
        #}}}
        #{{{ then deal with the units
        #}}}
        #{{{ finally, if I need to, reorder again to put the new dimension where I want it
        #}}}
        return self
    def chunk(self,axisin,*otherargs):
        r'''chunks axisin into multiple new axes
        arguments:
            axesout -- gives the names of the output axes
            shapesout -- optional -- if not given, it assumes equal length -- if given, one of the values can be -1, which is assumed length'''
        if len(otherargs) == 2:
            axesout,shapesout = otherargs
        elif len(otherargs) == 1:
            axesout = otherargs[0]
            shapesout = ndshape(self)[axisin]**(1./len(axesout))
            if abs(shapesout-round(shapesout)) > 1e-15: # there is some kind of roundoff error here
                raise ValueError('''In order for chunk to be called with
                        only a list of axes, the shape of the dimension you are
                        trying to split (here %s) must be an Nth root of
                        the original dimension size (here: %d), where N (here
                        %d) is the number of dimensions you are trying to chunk into'''%(axisin,ndshape(self)[axisin],len(axesout)))
            else:
                shapesout = round(shapesout)
            shapesout = [shapesout] * len(axesout)
        else:
            raise ValueError("otherargs must be one or two arguments!")
        if any([j == -1 for j in shapesout]):
            j = shapesout.index(-1)
            if j < len(shapesout)-1:
                shapesout[j] = int(round(ndshape(self)[axisin]/prod(r_[shapesout[0:j],shapesout[j+1:]])))
            else:
                shapesout[j] = int(round(ndshape(self)[axisin]/prod(shapesout[0:j])))
        if prod(shapesout) != ndshape(self)[axisin]:
            raise CustomError("The size of the axis (%s) you're trying to split (%s) doesn't match the size of the axes you're trying to split it into (%s = %s)"%(repr(axisin),repr(ndshape(self)[axisin]),repr(axesout),repr(shapesout)))
        thisaxis = self.axn(axisin)
        if self.getaxis(axisin) is not None and len(self.getaxis(axisin)) > 0:
            raise CustomError("You cannot chunk data with labels! Try 'chunk_by' instead!")
        #{{{ if there is a list of axis coordinates, add in slots for the new axes
        if type(self.axis_coords) is list:
            if len(self.axis_coords) == 0:
                self.axis_coords = [None] * len(self.dimlabels)
            for j in range(len(axesout)-1):
                self.axis_coords.insert(thisaxis,None)
        if type(self.axis_coords_error) is list:
            if len(self.axis_coords_error) == 0:
                self.axis_coords_error = [None] * len(self.dimlabels)
            for j in range(len(axesout)-1):
                self.axis_coords_error.insert(thisaxis,None)
        if type(self.axis_coords_units) is list:
            if len(self.axis_coords_units) == 0:
                self.axis_coords_units = [None] * len(self.dimlabels)
            for j in range(len(axesout)-1):
                self.axis_coords_units.insert(thisaxis,None)
        #}}}
        newshape = list(self.data.shape[0:thisaxis]) + shapesout + list(self.data.shape[thisaxis+1:])
        newnames = list(self.dimlabels[0:thisaxis]) + axesout + list(self.dimlabels[thisaxis+1:])
        self.data = self.data.reshape(newshape)
        self.dimlabels = newnames
        return self
    def chunkoff(self,axisin,newaxes,newshapes):
        r'''chunks up axisin, dividing it into newaxes with newshapes on the inside'''
        axesout = [axisin]+newaxes
        shapesout = [ndshape(self)[axisin]/prod(newshapes)]+newshapes
        return self.chunk(axisin,axesout,shapesout)
    def chunk_auto(self,axis_name,which_field,verbose = False,dimname = None):
        r'''assuming that axis "axis_name" is currently labeled with a structured array, choose one field ("which_field") of that structured array to generate a new dimension
        Note that for now, by definition, no error is allowed on the axes.
        However, once I upgrade to using structured arrays to handle axis and data errors, I will want to deal with that appropriately here.'''
        axis_number = self.axn(axis_name)
        new_axis,indices = unique(self.getaxis(axis_name)[which_field],
                return_inverse = True) # we are essentially creating a hash table for the axis.  According to numpy documentation, the hash indices that this returns should also be sorted sorted.
        if verbose: print "(chunk auto) indices look like this:",indices
        #{{{ check that there are equal numbers of all the unique new_axis
        index_count = array([count_nonzero(indices == j) for j in range(indices.max()+1)])
        if all(index_count == index_count[0]):
            if verbose: print "(chunk auto) Yes, there are equal numbers of all unique new_axis! (Each element of the hash table has been indexed the same number of times.)"
            #}}}
            #{{{ store the old shape and generate the new shape
            current_shape = list(self.data.shape)
            if verbose: print "(chunk auto) old shape -- ",current_shape
            new_shape = insert(current_shape,axis_number + 1,len(new_axis))
            new_shape[axis_number] /= len(new_axis) # the indices of the hash table become the new dimension
            #}}}
            #{{{ actually reorder the data and error -- perhaps a view would be more efficient here
            old_data = self.data
            has_data_error = not (self.get_error() is None)
            self.data = empty(new_shape,dtype = self.data.dtype)
            if has_data_error:
                old_error = self.get_error()
                self.set_error(empty(new_shape,dtype = self.data.dtype))
            #}}}
            #{{{ adjust all the relevant axis information
            #{{{ generate an axis label along the axis I'm chunking that's stripped of the field that I'm creating a dimension from (i.e. chunking off a new dimension based on) -- because I am independently manipulating the data, I don't use self.getaxis()
            x_strip_current_field = self.axis_coords[axis_number][
                    [j for j in
                        self.axis_coords[axis_number].dtype.names
                        if j != which_field]]
            #}}}
            #{{{ reshape the axis coordinate so that it becomes a 2D array with the new dimension chunked off
            self.axis_coords[axis_number] = empty((len(x_strip_current_field)/len(new_axis),len(new_axis)),
                    dtype = x_strip_current_field.dtype)
            if not (self.get_error(axis_name) is None):
                raise ValueError("Until I do the structured array upgrade chunk_auto will not be able to deal with an axis that has errors.")
            #}}}
            #{{{ everything is now ready to sort the data and residual axis into ordered slots
            #{{{ initialize the slices
            copy_to_slice   = len(new_shape)*[slice(None,None,None)] # this is the memory address inside the new data (where stuff goes)
            copy_from_slice = len(current_shape)*[slice(None,None,None)] # this is the memory address inside the old data (where stuff comes from)
            #}}}
            if has_data_error:
                data_error_location = self.get_error()
            for j in range(len(new_axis)): # j is the index in the hash table
                copy_to_slice[axis_number + 1]     = j
                copy_from_slice[axis_number]       = where(indices == j)[0]
                self.data[copy_to_slice]           = old_data[copy_from_slice]
                if has_data_error:
                    data_error_location[copy_to_slice] = old_error[copy_from_slice]
                if verbose: print "(chunk auto) ",j,'matches at',x_strip_current_field[copy_from_slice[axis_number]]
                self.axis_coords[axis_number][:,j] = x_strip_current_field[copy_from_slice[axis_number]]
            #}}}
            if verbose: print "(chunk auto) new axis -- ",self.axis_coords[axis_number]
            if verbose: print "(chunk auto) new shape -- ",self.data.shape
            #{{{ housekeeping for the various axes + data properties -- should perhaps be possible to do this first, then use .getaxis()
            self.dimlabels.insert(axis_number + 1,which_field)
            self.axis_coords.insert(axis_number + 1,new_axis)
            #{{{ by definition, axis can have neither errors nor units associated, for now.
            self.axis_coords_error.insert(axis_number + 1,None)
            self.axis_coords_units.insert(axis_number + 1,None)
            #}}}
            if verbose: print '(chunk auto) the dimensions of ',self.dimlabels[axis_number],'are (?? x ',self.dimlabels[axis_number+1],')=',self.axis_coords[axis_number].shape
            #}}}
            #}}}
            #{{{ deal appropriately with the "remainder axis" (axis_number)
            if dimname is None:
                remainder_axis_name = '_and_'.join(self.axis_coords[axis_number].dtype.names)
            else:
                remainder_axis_name = dimname
            #{{{ if everything is the same along the dimension that I've just
            # created (which is the second dimension), then get rid of the
            # duplicate labels
            test_axis = self.axis_coords[axis_number].T
            if verbose: print "(chunk auto) test axis -- ",test_axis
            test_axis = ascontiguousarray(test_axis).flatten().view([('',test_axis.dtype)]*test_axis.shape[1])
            if all(test_axis == test_axis[0]):
                self.axis_coords[axis_number] = self.axis_coords[axis_number][:,0].reshape(1,-1)
                if verbose: print "(chunk auto) collapsed to", self.axis_coords[axis_number]
            #}}}
            if self.axis_coords[axis_number].shape[0] == 1:# then this is a "valid" axis -- because, for each position of the new axis, there is only one value of the remainder axis
                self.axis_coords[axis_number] = self.axis_coords[axis_number].reshape(-1)
                self.dimlabels[axis_number] = remainder_axis_name
                if len(self.axis_coords[axis_number].dtype) == 1: # only one field, which by the previous line will be named appropriately, so drop the structured array name
                    new_dtype = self.axis_coords[axis_number].dtype.descr[0][1]
                    self.axis_coords[axis_number] = array(self.axis_coords[axis_number],dtype = new_dtype) # probably more efficiently done with a view, but leave alone for now
                return self
            else:
                #{{{ generate an index list to label the remainder axis, and generate a new nddata with the actual values (which are note copied across the new dimension that matches which_field)
                remainder_axis_index_list = r_[0:self.axis_coords[axis_number].shape[0]]
                new_data = nddata(self.axis_coords[axis_number],
                        self.axis_coords[axis_number].shape,
                        [remainder_axis_name,which_field])
                self.axis_coords[axis_number] = remainder_axis_index_list
                new_data.labels([remainder_axis_name,which_field],
                        [self.axis_coords[axis_number].copy(),self.axis_coords[axis_number + 1].copy()])
                self.dimlabels[axis_number] = remainder_axis_name + '_INDEX'
                #}}}
                return self,new_data
            #}}}
        else:
            raise ValueError("Along the axis '"+axis_name+"', the field '"+which_field+"' does not represent an axis that is repeated one or more times!  The counts for how many times each element along the field is used is "+repr(index_count))
            return
    #}}}
    #{{{ messing with data -- get, set, and copy
    def __getslice__(self,*args):
        print 'getslice! ',args
    def __setitem__(self,*args):
        righterrors = None
        A = args[0]
        if type(A) is nddata:
            _,B = self.aligndata(A)
            A = B.data # now the next part will handle this
        if type(A) is ndarray:# if selector is an ndarray
            if A.dtype is not dtype('bool'):
                raise ValueError("I don't know what to do with an ndarray subscript that has dtype "+repr(A.dtype))
            if A.shape != self.data.shape:
                temp = array(A.shape) == 1
                if all( array(A.shape)[temp] == array(self.data.shape)[temp]):
                    pass
                else:
                    raise ValueError("The shape of your logical mask "+repr(A.shape)+" and the shape of your data "+repr(self.data.shape)+" are not compatible (matching or singleton) -- I really don't think that you want to do this!")
            self.data[A] = args[1]
            return
        if isinstance(args[1],nddata):
            #{{{ reorder so the shapes match
            unshared_indices = list(set(args[1].dimlabels) ^ set(self.dimlabels))
            shared_indices = list(self.dimlabels)
            for j in unshared_indices:
                shared_indices.remove(j)
            if len(args[1].dimlabels) != len(shared_indices) or (not all([args[1].dimlabels[j] == shared_indices[j] for j in range(0,len(shared_indices))])):
                args[1].reorder[shared_indices]
            #}}}
            rightdata = args[1].data
            righterrors = args[1].get_error()
        else: # assume it's an ndarray
            rightdata = args[1]
            #{{{ if I just passed a function, assume that I'm applying some type of data-based mask
            if type(args[0]) is type(emptyfunction):
                thisfunc = args[0]
                self.data[thisfunc(self.data)] = rightdata
                return
                #}}}
            if (type(rightdata) is not ndarray): # in case its a scalar
                rightdata = array([rightdata])
        slicedict,axesdict,errordict,unitsdict = self._parse_slices(args[0]) # pull left index list from parse slices
        leftindex = self.fld(slicedict)
        rightdata = rightdata.squeeze()
        if len(rightdata.shape) > 0:
            left_shape = shape(self.data[tuple(leftindex)])
            try:
                self.data[tuple(leftindex)] = rightdata.reshape(left_shape) # assign the data
            except:
                raise CustomError('ERROR ASSIGNING NDDATA:\n',
                        'self.data.shape:',self.data.shape,
                        'left index',leftindex,'\n',
                        'rightdata.shape:',rightdata.shape,
                        '--> shape of left slice: ',left_shape)
        else:
            self.data[tuple(leftindex)] = rightdata
        lefterror = self.get_error()
        if lefterror is not None:
            lefterror[tuple(leftindex)] = righterrors.squeeze()
    def copy(self):
        return deepcopy(self)
    def __getitem__(self,args):
        if type(args) is type(emptyfunction):
            #{{{ just a lambda function operates on the data
            thisfunc = args
            newdata = self.copy()
            mask = thisfunc(newdata.data)
            newdata.data = newdata.data[mask]
            if len(newdata.dimlabels) == 1:
                x = newdata.getaxis(newdata.dimlabels[0])
                newdata.setaxis(newdata.dimlabels[0],x[mask])
            else:
                raise ValueError("I don't know how to do this for multidimensional data yet!")
            return newdata
            #}}}
        elif type(args) is nddata:
            #{{{ try the nddata mask
            A = args
            if isinstance(A,nddata) and A.data.dtype is dtype("bool"):
                thisshape = ndshape(A)
                nonsingleton = []
                for thisdim in A.dimlabels:
                    if thisshape[thisdim] != 1:
                        nonsingleton.append(thisdim)
                if len(nonsingleton) != 1:
                    raise ValueError("To index with an nddata, you must have only one dimension")
                else:
                    self.setaxis(nonsingleton[0],self.getaxis(nonsingleton[0])[A.data.flatten()])
                _,B = self.aligndata(A)
                A = B.data # now the next part will handle this
                if A.dtype is not dtype('bool'):
                    raise ValueError("I don't know what to do with an ndarray subscript that has dtype "+repr(A.dtype))
                if A.shape != self.data.shape:
                    temp = array(A.shape) == 1
                    if all( array(A.shape)[temp] == array(self.data.shape)[temp]):
                        pass
                    else:
                        raise ValueError("The shape of your logical mask "+repr(A.shape)+" and the shape of your data "+repr(self.data.shape)+" are not compatible (matching or singleton) -- I really don't think that you want to do this!")
                self.data = self.data[A]
                return self
            else:
                errmsg = "you passed a single argument of type "+repr(type(A))
                if type(A) is nddata:
                    errmsg += " with dtype "+repr(A.data.dtype)
                errmsg += " -- I don't know what to do with this"
                raise ValueError(errmsg)
            #}}}
        else:
            try:
                slicedict,axesdict,errordict,unitsdict = self._parse_slices(args)
            except:
                raise CustomError('error trying to get slices given by',args)
            if type(args) is not slice and type(args[1]) is list and type(args[0]) is str and len(args) == 2:
                return concat([self[args[0],x] for x in args[1]],args[0])
            indexlist = tuple(self.fld(slicedict))
            newlabels = [x for x in self.dimlabels if not isscalar(slicedict[x])] # generate the new list of labels, in order, for all dimensions that are not indexed by a scalar
        #{{{ properly index the data error
        if self.data_error != None:
            try:
                newerror = self.data_error[indexlist]
            except:
                raise ValueError('Problem trying to index data_error'+repr(self.data_error)+' with',repr(indexlist))
        else:
            newerror = None
        #}}}
        if len(self.axis_coords)>0:
            if errordict != None:
                axis_coords_error = [errordict[x] for x in newlabels]
            else:
                axis_coords_error = None
            if unitsdict is not None:
                axis_coords_units = [unitsdict[x] for x in newlabels]
            else:
                axis_coords_units = None
            try:
                retval =  nddata(self.data[indexlist],
                        self.data[indexlist].shape,
                        newlabels,
                        axis_coords = [axesdict[x] for x in newlabels],
                        axis_coords_error = axis_coords_error,
                        data_error = newerror,
                        other_info = self.other_info)
            except:
                raise CustomError("likely some problem recasting the data when trying to initialize a new nddata: shape of self.data",self.data,"indexlist",indexlist)
            retval.axis_coords_units = axis_coords_units
            retval.data_units = self.data_units
            return retval
        else:
            retval = nddata(self.data[indexlist],self.data[indexlist].shape,newlabels,self.other_info)
            retval.axis_coords_units = self.axis_coords_units
            retval.data_units = self.data_units
            return retval
    def _possibly_one_axis(self,*args):
        if len(args) == 1:
            return args[0]
        if len(args) > 1:
            raise ValueError('you can\'t pass more than one argument!!')
        if len(args) == 0:
            if len(self.dimlabels) == 1:
                axes = self.dimlabels
            elif len(self.dimlabels) == 0:
                raise ValueError("You're trying to do something to data with no dimensions")
            else:
                raise ValueError("If you have more than one dimension, you need to tell me which one!!")
        return axes
    def _process_kwargs(self,listoftuples,**kwargs):
        kwargnames,kwargdefaultvals = zip(*listoftuples)
        output = []
        for j,val in enumerate(kwargnames):
            output.append(kwargdefaultvals[j])
            if val in kwargs.keys():
                output[-1] = kwargs.pop(val)
        if len(kwargs) > 0:
            raise ValueError("I didn't understand the kwargs:",repr(kwargs))
        return tuple(output)
    def _parse_slices(self,args):
        """This controls nddata slicing:
            it previously took
            \'axisname\',value
            pairs where value was an index or a lambda function.
            Now, it also takes
            \'axisname\':value
            and
            \'axisname\':(value1,value2)
            pairs, where the values give either a single value or an inclusive range on the axis, respectively"""
        #print "DEBUG getitem called with",args
        errordict = None # in case it's not set
        if self.axis_coords_units is not None:
            unitsdict = self.mkd(self.axis_coords_units)
        axesdict = None # in case it's not set
        if type(args) is slice:
            args = [args]
        else:
            args = list(args)
        #{{{ to make things easier, convert "slice" arguments to "axis,slice" pairs
        j=0
        trueslice = [] # this lets me distinguish from 'axisname',slice
        while j < len(args):
            if type(args[j]) is slice:
                temp = args[j].start
                args.insert(j,temp)
                trueslice.append(temp)
                j+=1
            j+=2 # even only
        #}}}
        for j in range(0,len(args),2):
            if type(args[j]) is str_:
                args[j] = str(args[j]) # on upgrading + using on windows, this became necessary, for some reason I don't understand
        if type(args) in [float,int32,int,double]:
            raise CustomError('You tried to pass just a nddata[',type(args),']')
        if type(args[0]) is str:
            #{{{ create a slicedict and errordict to store the slices
            slicedict = dict(zip(list(self.dimlabels),[slice(None,None,None)]*len(self.dimlabels))) #initialize to all none
            if len(self.axis_coords)>0:
                #print "DEBUG --> trying to make dictionaries from axis coords of len",len(self.axis_coords),"and axis_coords_error of len",len(self.axis_coords_error),"when dimlabels has len",len(self.dimlabels)
                axesdict = self.mkd(self.axis_coords)
                errordict = self.mkd(self.axis_coords_error)
            for x,y in zip(args[0::2],args[1::2]):
                slicedict[x] = y
            #}}}
            #{{{ map the slices onto the axis coordinates and errors
            #print "DEBUG slicedict is",slicedict
            testf = lambda x: x+1
            if len(self.axis_coords)>0:
                for x,y in slicedict.iteritems():
                    #print "DEBUG, type of slice",x,"is",type(y)
                    if isscalar(y):
                        if axesdict[x] is not None:
                            axesdict.pop(x) # pop the axes for all scalar dimensions
                    elif type(y) is type(testf):
                        mask = y(axesdict[x])
                        slicedict[x] = mask
                        if axesdict[x] is not None:
                            axesdict[x] = axesdict[x][mask]
                    else:
                        if type(y) is slice and x in trueslice:
                            #print "DEBUG, I found",y,"to be of type slice"
                            if y.step is not None:
                                raise ValueError("setting the slice step is not currently supported")
                            else:
                                if type(y.stop) is tuple: #then I passed a single index
                                    temp = diff(axesdict[x]) 
                                    if not all(temp*sign(temp[0])>0):
                                        raise ValueError("you can only use the range format on data where the axis is in consecutively increasing or decreasing order")
                                    del temp
                                    if len(y.stop) > 2:
                                        raise ValueError("range with more than two values not currently supported")
                                    elif len(y.stop) == 1:
                                        temp_low = y.stop[0]
                                        temp_high = inf
                                        temp_high_value = inf
                                    else:
                                        temp_low = y.stop[0]
                                        temp_high = y.stop[1]
                                        temp_high_value = y.stop[1]
                                        if temp_low is None:
                                            temp_low = -inf
                                        if temp_high is None:
                                            temp_high = inf
                                        if temp_low > temp_high:
                                            temp_low,temp_high = temp_high,temp_low
                                    #print "DEBUG: slice values",temp_low,'to',temp_high
                                    if temp_low == inf:
                                        temp_low = axesdict[x].argmax()
                                    elif temp_low == -inf:
                                        temp_low = axesdict[x].argmin()
                                    else:
                                        temp_low = abs(axesdict[x] - temp_low).argmin()
                                    if temp_high == inf:
                                        temp_high = axesdict[x].argmax()
                                    elif temp_high == -inf:
                                        temp_high = axesdict[x].argmin()
                                    else:
                                        temp_high = abs(axesdict[x] - temp_high).argmin()
                                    if temp_high + 1 < len(axesdict[x]):
                                        #print "DEBUG: I test against value",axesdict[x][temp_high + 1]
                                        if axesdict[x][temp_high + 1] <= temp_high_value:
                                            temp_high += 1
                                    #print "DEBUG: evaluate to indices",temp_low,'to',temp_high,'out of',len(axesdict[x])
                                    if temp_high-temp_low == 0:
                                        raise ValueError('The indices '+repr(y)+' on axis '+x+' slices to nothing!  The limits of '+x+' are '+repr(axesdict[x].min())+':'+repr(axesdict[x].max()))
                                    slicedict[x] = slice(temp_low,temp_high,None)
                                    y = slicedict[x]
                                else: #then I passed a single index
                                    temp = abs(axesdict[x] - y.stop).argmin()
                                    slicedict[x] = slice(temp,temp+1,None)
                                    y = slicedict[x]
                        if axesdict[x] == []:
                            axesdict[x] = None
                        if axesdict[x] is not None:
                            try:
                                axesdict[x] = axesdict[x][y]
                            except:
                                raise ValueError("axesdict is"+repr(axesdict)+"and I want to set "+repr(x)+" subscript to its "+repr(y)+" value")
                if errordict != None and errordict != array(None):
                    for x,y in slicedict.iteritems():
                        if errordict[x] != None:
                            if isscalar(y):
                                errordict.pop(x)
                            elif type(y) is type(emptyfunction):
                                mask = y(axesdict[x])
                                errordict[x] = errordict[x][mask]
                            else:
                                try:
                                    errordict[x] = errordict[x][y] # default
                                except:
                                    raise CustomError('Trying to index',errordict,'-->',x,'=',errordict[x],'with',y,'error started as',self.axis_coords_error)
            if unitsdict != None and unitsdict != array(None):
                for x,y in slicedict.iteritems():
                    if unitsdict[x] != None:
                        if isscalar(y):
                            unitsdict.pop(x)
            return slicedict,axesdict,errordict,unitsdict
            #}}}
        else:
            raise CustomError('label your freaking dimensions! (type of args[0] is ',type(args[0]),'and it should be str!)')
    #}}}

    #{{{ hdf5 write
    def hdf5_write(self,h5path,verbose = False):
        #{{{ add the final node based on the name stored in the nddata structure
        if h5path[-1] != '/': h5path += '/' # make sure it ends in a slash first
        try:
            thisname = self.get_prop('name')
        except:
            raise CustomError("You're trying to save an nddata object which does not yet have a name, and you can't do this! Run yourobject.name('setname')")
        if type(thisname) is str:
            h5path += thisname
        else:
            raise CustomError("problem trying to store HDF5 file; you need to set the ``name'' property of the nddata object to a string first!")
        h5file,bottomnode = h5nodebypath(h5path) # open the file and move to the right node
        #print 'DEBUG 1: bottomnode is',bottomnode
        #}}}
        #{{{ print out the attributes of the data
        myattrs = normal_attrs(self)
        #{{{ separate them into data and axes
        mydataattrs = filter((lambda x: x[0:4] == 'data'),myattrs)
        myotherattrs = filter((lambda x: x[0:4] != 'data'),myattrs)
        myaxisattrs = filter((lambda x: x[0:4] == 'axis'),myotherattrs)
        myotherattrs = filter((lambda x: x[0:4] != 'axis'),myotherattrs)
        if verbose: print lsafe('data attributes:',zip(mydataattrs,map(lambda x: type(self.__getattribute__(x)),mydataattrs))),'\n\n'
        if verbose: print lsafe('axis attributes:',zip(myaxisattrs,map(lambda x: type(self.__getattribute__(x)),myaxisattrs))),'\n\n'
        if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
        #}}}
        #}}}
        #{{{ write the data table
        if 'data' in mydataattrs:
            if 'data_error' in mydataattrs and self.get_error() is not None and len(self.get_error()) > 0:
                thistable = rec.fromarrays([self.data,self.get_error()],names='data,error')
                mydataattrs.remove('data_error')
            else:
                thistable = rec.fromarrays([self.data],names='data')
            mydataattrs.remove('data')
            datatable = h5table(bottomnode,'data',thistable)
            #print 'DEBUG 2: bottomnode is',bottomnode
            #print 'DEBUG 2: datatable is',datatable
            if verbose: print "Writing remaining axis attributes\n\n"
            if len(mydataattrs) > 0:
                h5attachattributes(datatable,mydataattrs,self)
        else:
            raise CustomError("I can't find the data object when trying to save the HDF5 file!!")
        #}}}
        #{{{ write the axes tables
        if 'axis_coords' in myaxisattrs:
            if len(self.axis_coords) > 0:
                #{{{ create an 'axes' node
                axesnode = h5child(bottomnode, # current node
                        'axes', # the child
                        verbose = False,
                        create = True)
                #}}}
                for j,axisname in enumerate(self.dimlabels): # make a table for each different dimension
                    myaxisattrsforthisdim = dict([(x,self.__getattribute__(x)[j])
                        for x in list(myaxisattrs) if len(self.__getattribute__(x)) > 0]) # collect the attributes for this dimension and their values
                    if verbose: print lsafe('for axis',axisname,'myaxisattrsforthisdim=',myaxisattrsforthisdim)
                    if 'axis_coords' in myaxisattrsforthisdim.keys() and myaxisattrsforthisdim['axis_coords'] is not None:
                        if 'axis_coords_error' in myaxisattrsforthisdim.keys() and myaxisattrsforthisdim['axis_coords_error'] is not None and len(myaxisattrsforthisdim['axis_coords_error']) > 0: # this is needed to avoid all errors, though I guess I could use try/except
                            thistable = rec.fromarrays([myaxisattrsforthisdim['axis_coords'],myaxisattrsforthisdim['axis_coords_error']],names='data,error')
                            myaxisattrsforthisdim.pop('axis_coords_error')
                        else:
                            thistable = rec.fromarrays([myaxisattrsforthisdim['axis_coords']],names='data')
                        myaxisattrsforthisdim.pop('axis_coords')
                    datatable = h5table(axesnode,axisname,thistable)
                    #print 'DEBUG 3: axesnode is',axesnode
                    if verbose: print "Writing remaining axis attributes for",axisname,"\n\n"
                    if len(myaxisattrsforthisdim) > 0:
                        h5attachattributes(datatable,myaxisattrsforthisdim.keys(),myaxisattrsforthisdim.values())
        #}}}
        #{{{ Check the remaining attributes.
        if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
        if verbose: print "Writing remaining other attributes\n\n"
        if len(myotherattrs) > 0:
            #print 'DEBUG 4: bottomnode is',bottomnode
            test = repr(bottomnode) # somehow, this prevents it from claiming that the bottomnode is None --> some type of bug?
            try:
                h5attachattributes(bottomnode,
                    [j for j in myotherattrs if not self._contains_symbolic(j)],
                    self)
            except:
                raise CustomError('Problem trying to attach attributes',myotherattrs,'of self to node',bottomnode)
            warnlist = [j for j in myotherattrs if (not self._contains_symbolic(j)) and type(self.__getattribute__(j)) is dict]
            #{{{ to avoid pickling, test that none of the attributes I'm trying to write are dictionaries or lists
            if len(warnlist) > 0:
                print "WARNING!!, attributes",warnlist,"are dictionaries!"
            warnlist = [j for j in myotherattrs if (not self._contains_symbolic(j)) and type(self.__getattribute__(j)) is list]
            if len(warnlist) > 0:
                print "WARNING!!, attributes",warnlist,"are lists!"
            #}}}
            if verbose: print lsafe('other attributes:',zip(myotherattrs,map(lambda x: type(self.__getattribute__(x)),myotherattrs))),'\n\n'
        #}}}
        h5file.close()
    #}}}


#}}}
