class nddata_hdf5 (nddata):
    def __repr__(self):
        if hasattr(self,'_node_children'):
            return repr(self.datanode)
        else:
            return nddata.__repr__(self)
    def __del__(self):
        if hasattr(self,'_node_children'):
            self.h5file.close()
            del self.h5file
            del self.datanode
        return
    def __init__(self,pathstring):
        self.pathstring = pathstring
        try:
            self.h5file,self.datanode = h5nodebypath(pathstring,check_only = True)
        except:
            raise IndexError("I can't find the node"+pathstring)
        self._init_datanode(self.datanode)
    def _init_datanode(self,datanode,verbose = False,**kwargs):
        datadict = h5loaddict(datanode)
        #{{{ load the data, and pop it from datadict
        try:
            datarecordarray = datadict['data']['data'] # the table is called data, and the data of the table is called data
            mydata = datarecordarray['data']
        except:
            raise CustomError("I can't find the nddata.data")
        try:
            kwargs.update({'data_error':datarecordarray['error']})
        except:
            if verbose: print "No error found\n\n"
        datadict.pop('data')
        #}}}
        #{{{ be sure to load the dimlabels
        mydimlabels = datadict['dimlabels']
        if len(mydimlabels) == 1:
            if len(mydimlabels[0]) == 1:
                mydimlabels = list([mydimlabels[0][0]]) # for some reason, think I need to do this for length 1
        #}}}
        #{{{ load the axes and pop them from datadict
        datadict.pop('dimlabels')
        if 'axes' in datadict.keys():
            myaxiscoords = [None]*len(mydimlabels)
            myaxiscoordserror = [None]*len(mydimlabels)
            for axisname in datadict['axes'].keys():
                try:
                    axisnumber = mydimlabels.index(axisname)
                except AttributeError:
                    raise CustomError('mydimlabels is not in the right format!\nit looks like this:\n',mydimlabels,type(mydimlabels))
                recordarrayofaxis = datadict['axes'][axisname]['data']
                myaxiscoords[axisnumber] = recordarrayofaxis['data']
                if 'error' in recordarrayofaxis.dtype.names:
                    myaxiscoordserror[axisnumber] = recordarrayofaxis['error']
                datadict['axes'][axisname].pop('data')
                for k in datadict['axes'][axisname].keys():
                    print lsafen("Warning, attribute",k,"of axis table",axisname,"remains, but the code to load this is not yet supported")
                datadict['axes'].pop(axisname)
            kwargs.update({"axis_coords":myaxiscoords})
            kwargs.update({"axis_coords_error":myaxiscoordserror})
        elif len(mydimlabels)>1:
            raise CustomError("The current version uses the axis labels to figure out the shape of the data\nBecause you stored unlabeled data, I can\'t figure out the shape of the data!!")
            # the reshaping this refers to is done below
        #}}}
        nddata.__init__(self,
                mydata,
                mydata.shape,
                mydimlabels,
                **kwargs)
        #{{{ reshape multidimensional data to match the axes
        if len(mydimlabels)>1:
            det_shape = []
            for thisdimlabel in mydimlabels:
                try:
                    temp = self.getaxis(thisdimlabel)
                except:
                    temp = -1 # no axis is given
                if type(temp) is ndarray:
                    temp = len(temp)
                det_shape.append(temp)
            self.data = self.data.reshape(tuple([len(self.getaxis(x)) for x in mydimlabels]))
        #}}}
        for remainingattribute in datadict.keys():
            self.__setattr__(remainingattribute,datadict[remainingattribute])
        self.h5file.close()
        del self.h5file
        del self.datanode
        return
