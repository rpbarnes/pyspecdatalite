
    #{{{ poly. fit
    ## This also needs to be unclassed and lmfitted.
    def polyfit(self,axis,order=1,force_y_intercept = None):
        'return the coefficients and the fit --> later, should probably branch this off as a new type of fit class'
        x = self.getaxis(axis).copy().reshape(-1,1)
        #{{{ make a copy of self with the relevant dimension second to last (i.e. rows)
        formult = self.copy()
        neworder = list(formult.dimlabels)
        neworder.pop(neworder.index(axis))
        if len(neworder) > 1:
            neworder = neworder[:-1] + [axis] + neworder[-1]
        else:
            neworder = [axis] + neworder
        formult.reorder(neworder)
        #}}}
        y = formult.data
        #{{{ now solve Lx = y, where x is appropriate for our polynomial
        startingpower = 0
        if force_y_intercept != None:
            startingpower = 1
        L =  concatenate([x**j for j in range(startingpower,order+1)],axis=1) # note the totally AWESOME way in which this is done!
        #print 'fitting to matrix',L
        if force_y_intercept != None:
            y -= force_y_intercept
        c = dot(pinv(L),y)
        fity = dot(L,c)
        if force_y_intercept != None:
            #print "\n\nDEBUG: forcing from",fity[0],"to"
            fity += force_y_intercept
            #print "DEBUG: ",fity[0]
            c = c_[force_y_intercept,c]
        #}}}
        #{{{ rather than have to match up everything, just drop the fit data into formult, which should be the same size, shape, etc
        formult.data = fity
        formult.set_error(None)
        #}}}
        return c,formult
    #}}}

#{{{ fitdata
class fitdata(nddata):
    """
    This should be rewritten to make use of lmfit and scipy.linearize
    """
    def __init__(self,*args,**kwargs):
        #{{{ manual kwargs
        fit_axis = None
        if 'fit_axis' in kwargs.keys():
            fit_axis = kwargs.pop('fit_axis')
        #}}}
        if isinstance(args[0],nddata):
            #print "DEBUG trying to transfer",args[0].axis_coords_error
            myattrs = normal_attrs(args[0])
            for j in range(0,len(myattrs)):
                self.__setattr__(myattrs[j],args[0].__getattribute__(myattrs[j]))
            #nddata.__init__(self,
            #        args[0].data,
            #        args[0].data.shape,
            #        args[0].dimlabels,
            #        axis_coords = args[0].axis_coords,
            #        ft_start_time = args[0].ft_start_time,
            #        data_error = args[0].data_error,
            #        axis_coords_error = args[0].axis_coords_error,
            #        axis_coords_units = args[0].axis_coords_units,
            #        data_units = args[0].data_units,
            #        other_info = args[0].other_info,
            #        **kwargs)
        else:
            #self.__base_init(*args,**kwargs)
            nddata.__init__(self,*args,**kwargs)
        if fit_axis is None:
            if len(self.dimlabels) == 1:
                fit_axis = self.dimlabels[0]
            else:
                raise CustomError("I can't figure out the fit axis!")
        self.fit_axis = fit_axis
        #{{{ in the class, only store the forced values and indices they are set to
        self.set_to = None
        self.set_indices = None
        self.active_indices = None
        #}}}
        return
    def parameter_derivatives(self,xvals,set = None,set_to = None,verbose = False):
        r'return a matrix containing derivatives of the parameters, can set dict set, or keys set, vals set_to'
        if verbose: print 'parameter derivatives is called!'
        if iscomplex(self.data.flatten()[0]):
            print lsafen('Warning, taking only real part of fitting data!')
        if type(set) is dict:
            set_to = set.values()
            set = set.keys()
        solution_list = dict([(self.symbolic_dict[k],set_to[j])
            if k in set
            else (self.symbolic_dict[k],self.output(k))
            for j,k in enumerate(self.symbol_list)]) # load into the solution list
        number_of_i = len(xvals)
        parameters = self._active_symbols()
        mydiff_sym = [[]] * len(self.symbolic_vars)
        x = self.symbolic_x
        fprime = zeros([len(parameters),number_of_i])
        for j in range(0,len(parameters)):
            thisvar = self.symbolic_dict[parameters[j]]
            mydiff_sym[j] = sympy.diff(self.symbolic_func,thisvar)
            #print r'$\frac{\partial %s}{\partial %s}=%s$'%(self.function_name,repr(thisvar),sympy.latex(mydiff).replace('$','')),'\n\n'
            try:
                mydiff = mydiff_sym[j].subs(solution_list)
            except:
                raise CustomError('error trying to substitute',mydiff_sym[j],'with',solution_list)
            try:
                fprime[j,:] = array([complex(mydiff.subs(x,xvals[k])) for k in range(0,len(xvals))])
            except ValueError:
                raise CustomError('Trying to set index',j, 'shape(fprime)',shape(fprime), 'shape(xvals)',shape(xvals),'the thing I\'m trying to compute looks like this',[mydiff.subs(x,xvals[k]) for k in range(0,len(xvals))])
            except:
                raise CustomError('Trying to set index',j, 'shape(fprime)',shape(fprime), 'shape(xvals)',shape(xvals))
        return fprime
    def parameter_gradient(self,p,x,y,sigma):
        r'this gives the specific format wanted by leastsq'
        # for now, I'm going to assume that it's not using sigma, though this could be wrong
        # and I could need to scale everything by sigma in the same way as errfunc
        return self.parameter_derivatives(x,set = self.symbol_list,set_to = p).T
    def analytical_covariance(self):
        covarmatrix = zeros([len(self._active_symbols())]*2)
        #{{{ try this ppt suggestion --> his V is my fprime, but 
        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis))
        dirproductform = False
        if dirproductform:
            sigma = self.get_error()
            f1 = fprime.shape[0]
            f2 = fprime.shape[1]
            fprime1 = fprime.reshape(f1,1,f2) # j index
            fprime2 = fprime.reshape(1,f1,f2) # k index
            fprime_prod = fprime1 * fprime2
            fprime_prod = fprime_prod.reshape(-1,f2).T # direct product form
            try:
                covarmat = dot(pinv(fprime_prod),(sigma**2).reshape(-1,1))
            except ValueError:
                raise CustomError('shape of fprime_prod',shape(fprime_prod),'shape of inverse',shape(pinv(fprime_prod)),'shape of sigma',shape(sigma))
            covarmatrix = covarmat.reshape(f1,f1)
            for l in range(0,f1): 
                for m in range(0,f1): 
                    if l != m:
                        covarmatrix[l,m] /= 2
        else:
            sigma = self.get_error()
            #covarmatrix = dot(pinv(f),
            #        dot(diag(sigma**2),pinv(f.T)))
            J = matrix(fprime.T)
            #W = matrix(diag(1./sigma**2))
            #S = matrix(diag(sigma**2))
            #if hasattr(self,'data_covariance'):
            #    print "covariance data is present"
            S = matrix(self.get_covariance())
            Omegainv = S**-1
            #S = matrix(diag(sigma**2))
            #G = matrix(diag(1./sigma))
            #G = S**(-1/2) # analog of the above
            #covarmatrix = ((J.T * W * J)**-1) * J.T * W
            print 'a'
            minimizer = inv(J.T * Omegainv * J) * J.T * Omegainv
            covarmatrix = minimizer * S * minimizer.T
            #covarmatrix = array(covarmatrix * S * covarmatrix.T)
            #covarmatrix = array((J.T * G.T * G * J)**-1 * J.T * G.T * G * S * G.T * G * J * (J.T * G.T * G * J)**-1)
            #try:
            #    betapremult = (J.T * Omegainv * J)**-1 * J.T * Omegainv
            #except:
            #    print 'from sigma','\n\n',diag(sigma**2),'\n\n','from covarmatrix','\n\n',S,'\n\n'
            #    raise CustomError('problem generating estimator (word?)')
            #covarmatrix = array( betapremult * S * betapremult.T)
        #print "shape of fprime",shape(fprime),"shape of fprime_prod",shape(fprime_prod),'sigma = ',sigma,'covarmat=',covarmatrix,'\n'
        #}}}
        # note for this code, that it depends on above code I later moved to  parameter_derivatives
        #for j in range(0,shape(covarmatrix)[0]):
        #    for k in range(0,shape(covarmatrix)[0]):
        #        #mydiff_second = sympy.diff(mydiff_sym[j],self.symbolic_vars[k]).subs(solution_list)
        #        #fdprime = array([mydiff_second.subs(x,xvals[l])/sigma[l] for l in range(0,len(xvals))]) # only divide by sigma once, since there is only one f
        #        #try:
        #        temp = 1.0/(fprime[j,:] * fprime[k,:])
        #        mask = isinf(temp)
        #        covarmatrix[j,k] = mean(sigma[~mask]**2 * temp[~mask])# + 2. * mean(fminusE * fdprime)
        #        #except:
        #        #    raise CustomError('Problem multiplying covarmatrix', 'shape(fprime[j,:])',shape(fprime[j,:]), 'shape(fminusE)',shape(fminusE), 'shape(fdprime)',shape(fdprime))
        #        #if j != k:
        #        #    covarmatrix[j,k] *= 2
        return covarmatrix
    def gen_symbolic(self,function_name):
        r'''generates the symbolic representations the function'''
        self.function_name = function_name
        self.symbolic_vars = map(sympy.var,self.symbol_list)
        self.symbolic_x = sympy.var(self.fit_axis)
        #print lsafen('test symbol_list=',self.symbol_list)
        self.symbolic_dict = dict(zip(self.symbol_list,self.symbolic_vars))
        #print lsafen('test symbolic_vars=',self.symbolic_vars)
        #print lsafen('test symbolic_x=',self.symbolic_x)
        if hasattr(self,'fitfunc_raw_symb'):
            self.symbolic_func = self.fitfunc_raw_symb(self.symbolic_vars,self.symbolic_x)
        else:
            self.symbolic_func = self.fitfunc_raw(self.symbolic_vars,self.symbolic_x)
        self.function_string = sympy.latex(self.symbolic_func).replace('$','')
        self.function_string = r'$' + self.function_name + '=' + self.function_string + r'$'
        return self
    def copy(self): # for some reason, if I don't override this with the same thing, it doesn't override
        namelist = []
        vallist = []
        for j in dir(self):
            if self._contains_symbolic(j):
                namelist.append(j)
                vallist.append(self.__getattribute__(j))
                self.__delattr__(j)
        new = deepcopy(self)
        for j in range(0,len(namelist)):
            new.__setattr__(namelist[j],vallist[j])
        for j in range(0,len(namelist)):
            self.__setattr__(namelist[j],vallist[j])
        return new
    def gen_indices(self,set,set_to):
        r'''pass this set and set\_to parameters, and it will return:
        indices,values,mask
        indices --> gives the indices that are forced
        values --> the values they are forced to
        mask --> p[mask] are actually active in the fit'''
        if type(set) is not list:
            set = [set]
        if type(set_to) is not list:
            set_to = [set_to]
        if len(set) != len(set_to):
            raise CustomError('length of set=',set,'and set_to',set_to,'are not the same!')
        set_indices = map(self.symbol_list.index,set) # calculate indices once for efficiency
        active_mask = ones(len(self.symbol_list),dtype = bool)
        active_mask[set_indices] = False # generate the mask of indices that are actively fit
        return set_indices,set_to,active_mask
    def remove_inactive_p(self,p):
        return p[self.active_mask]
    def add_inactive_p(self,p):
        if self.set_indices != None:
            #{{{ uncollapse the function
            temp = p.copy()
            p = zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indices] = self.set_to # then just set the forced values to their given values
        return p
    def fitfunc(self,p,x):
        r"this wraps fitfunc_raw (which gives the actual form of the fit function) to take care of forced variables"
        p = self.add_inactive_p(p)
        return self.fitfunc_raw(p,x)
    def errfunc(self,p,x,y,sigma):
        '''just the error function'''
        fit = self.fitfunc(p,x)
        #normalization = sum(1.0/sigma)
        #print 'DEBUG: y=',y,'\nfit=',fit,'\nsigma=',sigma,'\n\n'
        sigma[sigma == 0.0] = 1
        try:
            retval = (y-fit)/sigma #* normalization
            #print 'DEBUG: retval=',retval,'\n\n'
        except ValueError:
            raise CustomError('your error (',shape(sigma),') probably doesn\'t match y (',shape(y),') and fit (',shape(fit),')')
        return retval
    def pinv(self,*args,**kwargs):
        if 'verbose' in kwargs.keys():
            verbose = kwargs.pop('verbose')
        else:
            verbose = False
        retval = self.linear(*args,**kwargs)
        y = retval.data
        yerr = retval.get_error()
        x_axis = retval.dimlabels[0]
        x = retval.getaxis(x_axis)
        nopowerindex = argmax(x)
        mask = logical_not(r_[0:len(x)] == nopowerindex)
        y = y[mask]
        yerr = yerr[mask]
        x = x[mask]
        L = c_[x.reshape((-1,1)),ones((len(x),1))]
        retval = dot(pinv(L,rcond = 1e-17),y)
        if verbose:
            print r'\label{fig:pinv_figure_text}y=',y,'yerr=',yerr,'%s='%x_axis,x,'L=',L
            print '\n\n'
            print 'recalc y = ',dot(L,retval)
            print 'recalc E = ',1.0-1.0/dot(L,retval)
            print 'actual E = ',self.data
        return retval
    def linear(self,*args,**kwargs):
        r'''return the linear-form function, either smoothly along the fit function, or on the raw data, depending on whether or not the taxis argument is given
        can take optional arguments and pass them on to eval'''
        #print "DEBUG called linear"
        if len(args) == 1:
            taxis = self._taxis(args[0]) # handle integer as well
            return self.linfunc(taxis,self.eval(taxis,**kwargs).data) # if we pass an argument, return the function across the entire time axis passed
        else:
            return self.linfunc(self.getaxis(self.fit_axis),self.data,yerr = self.get_error(),xerr = self.get_error(self.fit_axis)) # otherwise, return the raw data
    def output(self,*name):
        r'''give the fit value of a particular symbol'''
        if not hasattr(self,'fit_coeff') or self.fit_coeff is None:
            return None
        p = self.fit_coeff.copy()
        if self.set_indices != None:
            #{{{ uncollapse the function
            temp = p.copy()
            p = zeros(len(self.symbol_list))
            p[self.active_mask] = temp
            #}}}
            p[self.set_indices] = self.set_to # then just set the forced values to their given values
            #print "DEBUG trying to uncollapse in fitfunc w/ ",self.symbol_list,"; from",temp,"to",p
        # this should also be generic
        if len(name) == 1:
            try:
                return p[self.symbol_list.index(name[0])]
            except:
                raise CustomError("While running output: couldn't find",name,"in",self.symbol_list)
        elif len(name) == 0:
            # return a record array
            return array(tuple(p),{"names":list(self.symbol_list),"formats":['double']*len(p)}).reshape(1)
        else:
            raise CustomError("You can't pass",len(name),"arguments to .output()")
    def _pn(self,name):
        return self.symbol_list.index(name)
    def _active_symbols(self):
        if not hasattr(self,'active_symbols'):
            if self.set_indices != None:
                self.active_symbols = [x for x in self.symbol_list if self.active_mask[self._pn(x)]]
            else:
                self.active_symbols = list(self.symbol_list)
        return self.active_symbols
    def _pn_active(self,name):
        return self._active_symbols().index(name)
    def covar(self,*names):
        r'''give the covariance for the different symbols'''
        if len(names) == 1:
            names = [names[0],names[0]]
        if self.covariance is not None:
            return self.covariance[self._pn_active(names[0]),
                    self._pn_active(names[1])].copy()
        else:
            return None
    def covarmat(self,*names):
        if (len(names) == 1) and (names[0] is 'recarray'):
            if hasattr(self,'active_mask'):
                active_symbols = [x for x in self.symbol_list if self.active_mask[self._pn(x)]]
            else:
                active_symbols = list(self.symbol_list)
            if len(active_symbols) != self.covariance.shape[0]:
                raise CustomError('length of active symbols',active_symbols,'doesnt match covariance matrix size(',self.covariance.shape[0],')!')
            recnames = ['labels'] + active_symbols
            recdata = []
            for j in range(0,self.covariance.shape[0]): 
                thisdata = [active_symbols[j]] + list(double(self.covariance[j,:].copy())) # the first index is the row
                recdata.append(make_rec(thisdata,recnames))
            return r_[tuple(recdata)]
        if len(names) > 0:
            indices = map(self._pn_active,names) # slice out only these rows and columns
            return self.covariance[r_[indices],:][:,r_[indices]].copy()
        else:
            try:
                return self.covariance.copy()
            except:
                return zeros([len(self.fit_coeff)]*2,dtype = 'double')
    def latex(self,verbose = False):
        r'''show the latex string for the function, with all the symbols substituted by their values'''
        # this should actually be generic to fitdata
        p = self.fit_coeff
        printfstring = self.function_string
        printfargs = []
        allsymb = []
        locations = []
        for j in range(0,len(self.symbol_list)):
            #symbol = self.symbol_list[j]
            symbol = sympy.latex(self.symbolic_vars[j]).replace('$','')
            if verbose: print 'DEBUG: replacing symbol \\verb|',symbol,'|'
            location = printfstring.find(symbol)
            while location != -1:
                if printfstring[location-1] == '-':
                    newstring = printfstring[:location-1]+'+%01.03g'+printfstring[location+len(symbol):] # replace the symbol in the written function with the appropriate number
                    thissign = -1.0
                else:
                    newstring = printfstring[:location]+'%01.03g'+printfstring[location+len(symbol):] # replace the symbol in the written function with the appropriate number
                    thissign = 1.0
                if verbose: print r"\begin{verbatim} trying to replace",printfstring[location:location+len(symbol)],r'\end{verbatim}'
                printfstring = newstring
                printfargs += [thissign*p[j]] # add that number to the printf list
                locations += [location]
                allsymb += [symbol]
                location = printfstring.find(symbol)
        printfargs = [printfargs[x] for x in argsort(locations)]
        if verbose: print r"\begin{verbatim}trying to generate",self.function_string,'\n',printfstring,'\n',[allsymb[x] for x in argsort(locations)],'\n',printfargs,r'\end{verbatim}'
        return printfstring%tuple(printfargs)
    def settoguess(self):
        'a debugging function, to easily plot the initial guess'
        self.fit_coeff = real(self.guess())
        return self
    def _taxis(self,taxis):
        r'You can enter None, to get the fit along the same range as the data, an integer to give the number of points, or a range of data, which will return with 300 points'
        if taxis is None:
            taxis = self.getaxis(self.fit_axis).copy()
        elif type(taxis) is int:
            taxis = linspace(self.getaxis(self.fit_axis).min(),
                    self.getaxis(self.fit_axis).max(),
                    taxis)
        elif not isscalar(taxis) and len(taxis) == 2:
            taxis = linspace(taxis[0],taxis[1],300)
        return taxis
    def eval(self,taxis,set = None,set_to = None):
        r'''after we have fit, evaluate the fit function along the axis taxis
        set and set_to allow you to forcibly set a specific symbol to a specific value --> however, this does not affect the class, but only the return value'''
        if type(set) is dict:
            set_to = set.values()
            set = set.keys()
        taxis = self._taxis(taxis)
        if hasattr(self,'fit_coeff') and self.fit_coeff is not None:
            p = self.fit_coeff.copy()
        else:
            p = array([NaN]*len(self.symbol_list))
        #{{{ LOCALLY apply any forced values
        if set != None:
            if self.set_indices != None:
                raise CustomError("your'e trying to set indices in an eval function for a function that was fit constrained; this is not currently supported")
            set_indices,set_to,active_mask = self.gen_indices(set,set_to)
            p[set_indices] = set_to
        #}}}
        #{{{ make a new, blank array with the fit axis expanded to fit taxis
        newdata = ndshape(self)
        newdata[self.fit_axis] = size(taxis)
        newdata = newdata.alloc()
        newdata.set_plot_color(self.get_plot_color())
        #}}}
        #{{{ keep all axis labels the same, except the expanded one
        newdata.axis_coords = list(newdata.axis_coords)
        newdata.labels([self.fit_axis],list([taxis]))
        #}}}
        newdata.data[:] = self.fitfunc(p,taxis).flatten()
        return newdata
    def makereal(self):
        self.data = real(self.data)
        return
    def rename(self,previous,new):
        if previous == self.fit_axis:
            self.fit_axis = new
        nddata.rename(self,previous,new)
        return self
    def fit(self,set = None, set_to = None, force_analytical = False):
        r'''actually run the fit'''
        if type(set) is dict:
            set_to = set.values()
            set = set.keys()
        x = self.getaxis(self.fit_axis)
        if iscomplex(self.data.flatten()[0]):
            print lsafen('Warning, taking only real part of fitting data!')
        y = real(self.data)
        sigma = self.get_error()
        if sigma is None:
            print '{\\bf Warning:} You have no error associated with your plot, and I want to flag this for now\n\n'
            warnings.warn('You have no error associated with your plot, and I want to flag this for now',Warning)
            sigma = ones(shape(y))
        p_ini = real(array(self.guess())) # need the numpy format to allow boolean mask
        if set != None:
            self.set_indices,self.set_to,self.active_mask = self.gen_indices(set,set_to)
            p_ini = self.remove_inactive_p(p_ini)
        leastsq_args = (self.errfunc, p_ini)
        leastsq_kwargs = {'args':(x,y,sigma),
                    'full_output':True}# 'maxfev':1000*(len(p_ini)+1)}
        if hasattr(self,'has_grad') and self.has_grad == True:
            leastsq_kwargs.update({'Dfun':self.parameter_gradient})
        if 'Dfun' in leastsq_kwargs.keys():
            print "yes, Dfun passed with arg",leastsq_kwargs['Dfun']
        try:
            p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
        #{{{ just give various explicit errors
        except TypeError,err:
            if type(x) != ndarray and type(y) != ndarray:
                raise CustomError('leastsq failed because the two arrays aren\'t of the right type','type(x):',type(x),'type(y):',type(y))
            else:
                if any(shape(x) != shape(y)):
                    raise CustomError('leastsq failed because the two arrays do not match in size size','shape(x):',shape(x),'shape(y):',shape(y))
            raise CustomError('leastsq failed because of a type error!','type(x):',showtype(x),'type(y):',showtype(y),'type(sigma)',showtype(sigma),'shape(x):',shape(x),'shape(y):',shape(y),'shape(sigma)',shape(sigma),'p_ini',type(p_ini),p_ini)
        except ValueError,err:
            raise CustomError('leastsq failed with "',err,'", maybe there is something wrong with the input:',self)
        except:
            raise CustomError('leastsq failed; I don\'t know why')
        #}}}
        if success not in [1,2,3,4]:
            #{{{ up maximum number of evals
            if mesg.find('maxfev'):
                leastsq_kwargs.update({ 'maxfev':50000 })
                p_out,cov,infodict,mesg,success = leastsq(*leastsq_args,**leastsq_kwargs)
                if success != 1:
                    if mesg.find('two consecutive iterates'):
                        print r'{\Large\color{red}{\bf Warning data is not fit!!! output shown for debug purposes only!}}','\n\n'
                        print r'{\color{red}{\bf Original message:}',lsafe(mesg),'}','\n\n'
                        infodict_keys = infodict.keys()
                        infodict_vals = infodict.values()
                        if 'nfev' in infodict_keys:
                            infodict_keys[infodict_keys.index('nfev')] = 'nfev, number of function calls'
                        if 'fvec' in infodict_keys:
                            infodict_keys[infodict_keys.index('fvec')] = 'fvec, the function evaluated at the output'
                        if 'fjac' in infodict_keys:
                            infodict_keys[infodict_keys.index('fjac')] = 'fjac, A permutation of the R matrix of a QR factorization of the final approximate Jacobian matrix, stored column wise. Together with ipvt, the covariance of the estimate can be approximated.'
                        if 'ipvt' in infodict_keys:
                            infodict_keys[infodict_keys.index('ipvt')] = 'ipvt, an integer array of length N which defines a permutation matrix, p, such that fjac*p = q*r, where r is upper triangular with diagonal elements of nonincreasing magnitude.  Column j of p is column ipvt(j) of the identity matrix'
                        if 'qtf' in infodict_keys:
                            infodict_keys[infodict_keys.index('qtf')] = 'qtf, the vector (transpose(q)*fvec)'
                        for k,v in zip(infodict_keys,infodict_vals):
                            print r'{\color{red}{\bf %s:}%s}'%(k,v),'\n\n'
                        #self.fit_coeff = None
                        #self.settoguess()
                        #return
                    else:
                        raise CustomError('leastsq finished with an error message:',mesg)
                    #}}}
            else:
                raise CustomError('leastsq finished with an error message:',mesg)
        else:
            print r'{\color{blue}'
            print lsafen("Fit finished successfully with a code of %d and a message ``%s''"%(success,mesg))
            print r'}'
        self.fit_coeff = p_out # note that this is stored in HIDDEN form
        dof = len(x) - len(p_out)
        if hasattr(self,'symbolic_x') and force_analytical:
            self.covariance = self.analytical_covariance()
        else:
            if force_analytical: raise CustomError("I can't take the analytical covariance!  This is problematic.")
            if cov == None:
                #raise CustomError('cov is none! why?!, x=',x,'y=',y,'sigma=',sigma,'p_out=',p_out,'success=',success,'output:',p_out,cov,infodict,mesg,success)
                print r'{\color{red}'+lsafen('cov is none! why?!, x=',x,'y=',y,'sigma=',sigma,'p_out=',p_out,'success=',success,'output:',p_out,cov,infodict,mesg,success),'}\n'
            self.covariance = cov
        if self.covariance is not None:
            try:
                self.covariance *= sum(infodict["fvec"]**2)/dof # scale by chi_v "RMS of residuals"
            except:
                raise CustomError("type(self.covariance)",type(self.covariance),
                        "type(infodict[fvec])",type(infodict["fvec"]),
                        "type(dof)",type(dof))
        #print lsafen("DEBUG: at end of fit covariance is shape",shape(self.covariance),"fit coeff shape",shape(self.fit_coeff))
        return
    def bootstrap(self,points,swap_out = exp(-1.0),seedval = 10347,minbounds = {},maxbounds = {}):
        print r'\begin{verbatim}'
        seed(seedval)
        fitparameters = list(self.symbol_list)
        recordlist = array([tuple([0]*len(fitparameters))]*points,
                {'names':tuple(fitparameters),'formats':tuple(['double']*len(fitparameters))}) # make an instance of the recordlist
        for runno in range(0,points):
            success = False # because sometimes this doesn't work
            while success is False:
                thiscopy = self.copy()
                #{{{ discard datapoints
                origsizecheck = double(size(thiscopy.data))
                mask = thiscopy.random_mask(thiscopy.fit_axis,threshold = swap_out)
                thiscopy.data = thiscopy.data[mask]
                derr = thiscopy.get_error()
                x = thiscopy.getaxis(thiscopy.fit_axis)
                x = x[mask] # note that x is probably no longer a pointer
                derr = derr[mask]
                #print 'DEBUG: size of data after cut',double(size(thiscopy.data))/origsizecheck,' (expected ',1.-swap_out,')'
                #}}}
                #{{{ now extend
                number_to_replace = origsizecheck - thiscopy.data.size
                #print 'DEBUG: number_to_replace',number_to_replace
                random_indices = int32((rand(number_to_replace)*(thiscopy.data.size-1.0)).round())
                thiscopy.data = r_[thiscopy.data,thiscopy.data.copy()[random_indices]]
                thiscopy.labels([thiscopy.fit_axis],[r_[x,x.copy()[random_indices]]])
                thiscopy.set_error(r_[derr,derr.copy()[random_indices]])
                #print 'DEBUG: size of data after extension',double(size(thiscopy.data))/origsizecheck
                #}}}
                try:
                    thiscopy.fit()
                    success = True
                    if len(minbounds) > 0:
                        for k,v in minbounds.iteritems():
                            if thiscopy.output(k) < v:
                                success = False
                    if len(maxbounds) > 0:
                        for k,v in maxbounds.iteritems():
                            if thiscopy.output(k) > v:
                                success = False
                except:
                    #print 'WARNING, didn\'t fit'
                    success = False
                # here, use the internal routines, in case there are constraints, etc
                if success is True:
                    for name in thiscopy.symbol_list: # loop over all fit coeff
                        recordlist[runno][name] = thiscopy.output(name)
        print r'\end{verbatim}'
        return recordlist # collect into a single recordlist array
    def guess(self,verbose = False,super_verbose = False):
        r'''provide the guess for our parameters; by default, based on pseudoinverse'''
        self.has_grad = False
        if iscomplex(self.data.flatten()[0]):
            print lsafen('Warning, taking only real part of fitting data!')
        y = real(self.data)
        # I ended up doing the following, because as it turns out
        # T1 is a bad fit function, because it has a singularity!
        # this is probably why it freaks out if I set this to zero
        # on the other hand, setting a value of one seems to be
        # bad for very short T1 samples
        which_starting_guess = 0
        thisguess = self.starting_guesses[which_starting_guess]
        numguesssteps = 20
        #{{{ for some reason (not sure) adding a dimension to y
        new_y_shape = list(y.shape)
        new_y_shape.append(1)
        y = y.reshape(tuple(new_y_shape))
        #}}}
        #{{{ evaluate f, fprime and residuals
        guess_dict = dict(zip(self.symbol_list,list(thisguess)))
        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
        f_at_guess = real(self.eval(None,set = guess_dict).data)
        try:
            f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
        except:
            raise CustomError('trying to reshape f_at_ini_guess from',f_at_guess.shape,'to',new_y_shape)
        thisresidual = sqrt((y-f_at_guess)**2).sum()
        #}}}
        lastresidual = thisresidual
        for j in range(0,numguesssteps):
            if super_verbose: print '\n\n(matlablike.guess) '+r'\begin{verbatim} fprime = \n',fprime,'\nf_at_guess\n',f_at_guess,'y=\n',y,'\n',r'\end{verbatim}'
            if super_verbose: print '\n\n(matlablike.guess) shape of parameter derivatives',shape(fprime),'shape of output',shape(y),'\n\n'
            regularization_bad = True
            alpha_max = 100.
            alpha_mult = 2.
            alpha = 0.1 # maybe I can rather estimate this based on the change in the residual, similar to in L-M?
            if verbose: print '\n\n(matlablike.guess) value of residual before regularization %d:'%j,thisresidual
            while regularization_bad:
                newguess = real(array(thisguess) + dot(pinvr(fprime.T,alpha),(y-f_at_guess)).flatten())
                mask = newguess < self.guess_lb
                newguess[mask] = self.guess_lb[mask]
                mask = newguess > self.guess_ub
                newguess[mask] = self.guess_ub[mask]
                if any(isnan(newguess)):
                    if verbose: print '\n\n(matlablike.guess) Regularization blows up $\\rightarrow$ increasing $\\alpha$ to %0.1f\n\n'%alpha
                    alpha *= alpha_mult
                else:
                    #{{{ evaluate f, fprime and residuals
                    guess_dict = dict(zip(self.symbol_list,list(newguess)))
                    # only evaluate fprime once we know this is good, below
                    f_at_guess = real(self.eval(None,set = guess_dict).data)
                    try:
                        f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                    except:
                        raise CustomError('trying to reshape f_at_ini_guess from',f_at_guess.shape,'to',new_y_shape)
                    thisresidual = sqrt((y-f_at_guess)**2).sum()
                    #}}}
                    if (thisresidual-lastresidual)/lastresidual > 0.10:
                        alpha *= alpha_mult
                        if verbose: print '\n\n(matlablike.guess) Regularized Pinv gave a step uphill $\\rightarrow$ increasing $\\alpha$ to %0.1f\n\n'%alpha
                    else: # accept the step
                        regularization_bad = False
                        thisguess = newguess
                        lastresidual = thisresidual
                        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
                if alpha > alpha_max:
                    print "\n\n(matlablike.guess) I can't find a new guess without increasing the alpha beyond %d\n\n"%alpha_max
                    if which_starting_guess >= len(self.starting_guesses)-1:
                        print "\n\n(matlablike.guess) {\\color{red} Warning!!!} ran out of guesses!!!%d\n\n"%alpha_max
                        return thisguess
                    else:
                        which_starting_guess += 1
                        thisguess = self.starting_guesses[which_starting_guess]
                        print "\n\n(matlablike.guess) try a new starting guess:",lsafen(thisguess)
                        j = 0 # restart the loop
                        #{{{ evaluate f, fprime and residuals for the new starting guess
                        guess_dict = dict(zip(self.symbol_list,list(thisguess)))
                        fprime = self.parameter_derivatives(self.getaxis(self.fit_axis),set = guess_dict)
                        f_at_guess = real(self.eval(None,set = guess_dict).data)
                        try:
                            f_at_guess = f_at_guess.reshape(tuple(new_y_shape))
                        except:
                            raise CustomError('trying to reshape f_at_ini_guess from',f_at_guess.shape,'to',new_y_shape)
                        thisresidual = sqrt((y-f_at_guess)**2).sum()
                        #}}}
                        regularization_bad = False # jump out of this loop
            if verbose: print '\n\n(matlablike.guess) new value of guess after regularization:',lsafen(newguess)
            if verbose: print '\n\n(matlablike.guess) value of residual after regularization:',thisresidual
        return thisguess
#}}}
