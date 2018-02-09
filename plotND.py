
    #{{{ 3D mesh plot All of the following need to be un-Classified
    def meshplot(self,stride = None,alpha = 0.3,onlycolor = False,light = None,rotation = None,cmap = cm.gray,ax = None,invert = False,**kwargs):
        r'''takes both rotation and light as elevation, azimuth
        only use the light kwarg to generate a black and white shading display'''
        X,Y,Z = self.matrices_3d()
        if light == True:
            light = [0,0]# I think this is 45 degrees up shining down from the left of the y axis
        if not onlycolor:
            if ax is None: 
                ax = self._init_3d_axis(ax,rotation = rotation)
            else:
                if rotation is not None:
                    raise ValueError("you can only set the rotation once! (you tried"+repr(rotation)+")")
        rstride = 1
        cstride = 1
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        if stride is not None:
            if x_dim in stride.keys():
                rstride = stride[x_dim]
            if y_dim in stride.keys():
                cstride = stride[y_dim]
        if light is not None:
            ls = LightSource(azdeg = light[1],altdeg = light[0])
            if cmap is not None:
                rgb = ls.shade(Z,cmap)
        else:
            mask = isfinite(Z.flatten())
            for_rgb = Z-Z.flatten()[mask].min()
            for_rgb /= for_rgb.flatten()[mask].max()
            if cmap is not None:
                rgb = cmap(for_rgb)
        if onlycolor:
            imshow(rgb)
        else:
            if light is None:
                if cmap is not None:
                    kwargs.update({'cmap':cmap})
                ax.plot_surface(X,Y,Z,
                        rstride = rstride,
                        cstride = cstride,
                        shade = True,
                        **kwargs)
            else:
                newkwargs = {}
                newkwargs['linewidth'] = 0.0
                newkwargs.update(kwargs)
                if cmap is not None:
                    newkwargs['facecolors'] = rgb
                ax.plot_surface(X,Y,Z,
                        rstride = rstride,
                        cstride = cstride,
                        alpha = alpha,
                        shade = False,
                        **newkwargs)
            ax.set_xlabel(x_dim)
            ax.set_ylabel(y_dim)
            ax.set_zlabel(self.name())
        if onlycolor:
            return
        else:
            return ax

    def contour(self,levels = True):
        x_axis,y_axis = self.dimlabels
        x = self.getaxis(x_axis)[:,None]
        y = self.getaxis(y_axis)[None,:]
        cs = contour(x*ones_like(y),ones_like(x)*y,self.data,levels = r_[self.data.min():self.data.max():30j])
        if levels:
            clabel(cs,inline = 1,fontsize = 10)
        xlabel(unitify_axis(self,x_axis))
        ylabel(unitify_axis(self,y_axis))
        return cs

    def waterfall(self,alpha = 0.3,ax = None,rotation = None,color = 'b',edgecolor = 'k'):
        if ax is None: 
            ax = self._init_3d_axis(ax,rotation = rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        if len(self.dimlabels) > 2:
            raise CustomError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = self.dimlabels[0]
        y_dim = self.dimlabels[1]
        try:
            x_axis = self.retaxis(x_dim).data
        except:
            raise CustomError('trying to get the info on axis',x_dim,'which is',self.getaxis(x_dim))
        y_axis = self.retaxis(y_dim).data
        #}}}
        ax.set_xlabel(unitify_axis(self,x_dim))
        ax.set_ylabel(unitify_axis(self,y_dim))
        ax.set_zlabel(unitify_axis(self,self.name(),is_axis = False))
        verts = []
        xs = x_axis.flatten()
        xs = r_[xs[0],xs,xs[-1]] # add points for the bottoms of the vertices
        ys = y_axis.flatten()
        for j in range(0,len(ys)):
            zs = self[y_dim,j].data.flatten()
            zs = r_[0,zs,0]
            verts.append(zip(xs,zs)) # one of the faces
        poly = PolyCollection(verts, facecolors = [color]*len(verts), edgecolors = edgecolor) # the individual facecolors would go here
        poly.set_alpha(alpha)
        fig = gcf()
        ax.add_collection3d(poly,zs = ys, zdir = 'y')
        ax.set_zlim3d(self.data.min(),self.data.max())
        ax.set_xlim3d(xs.min(),xs.max())
        ax.set_ylim3d(ys.min(),ys.max())
        return ax

    def _init_3d_axis(self,ax,rotation = None):
        # other things that should work don't work correctly, so use this to initialize the 3D axis
        #ax.view_init(elev = rotation[0],azim = rotation[1])
        if rotation is None:
            rotation = [0,0]
        if ax == None:
            fig = gcf()
            ax = axes3d.Axes3D(fig)
            print "I'm trying to rotate to",rotation
            #ax.view_init(20,-120)
            #ax.view_init(elev = 20 + rotation[1],azim = -120 + rotation[0])
            ax.view_init(azim = rotation[0],elev = rotation[1])
        return ax

    def oldtimey(self,alpha = 0.5,ax = None,linewidth = None,sclinewidth = 20.,light = True,rotation = None,invert = False,**kwargs):
        sortedself = self.copy()
        self.sort(self.dimlabels[0])
        self.sort(self.dimlabels[1])
        if invert:
            print "trying to invert oldtimey"
        if linewidth == None:
            linewidth = sclinewidth/sortedself.data.shape[1]
            print "setting linewidth to %0.1f"%linewidth
        if ax is None: 
            ax = sortedself._init_3d_axis(ax,rotation = rotation)
        else:
            if rotation is not None:
                raise ValueError("you can only set the rotation once!")
        ax = sortedself.meshplot(linewidth = 0,light = light,ax = ax,invert = invert)
        #return
        if len(sortedself.dimlabels) > 2:
            raise CustomError("I don't know how to handle something with more than two dimensions for a surface plot!")
        #{{{ shared to both
        x_dim = sortedself.dimlabels[0]
        y_dim = sortedself.dimlabels[1]
        x_axis = sortedself.retaxis(x_dim).data
        y_axis = sortedself.retaxis(y_dim).data
        #}}}
        verts = []
        xs = x_axis.flatten()
        ys = y_axis.flatten() # this is the depth dimension
        if invert:
            ys = ys[::-1]
        for j in range(0,len(ys)):
            zs = sortedself[y_dim,j].data.flatten() # pulls the data (zs) for a specific y slice
            if invert:
                zs = zs[::-1]
            ax.plot(xs,ones(len(xs))*ys[j],zs,'k',linewidth = linewidth)
        fig = gcf()
        ax.set_zlim3d(sortedself.data.min(),sortedself.data.max())
        ax.set_xlim3d(xs.min(),xs.max())
        #if invert:
        #    ax.set_ylim3d(ys.max(),ys.min())
        #else:
        ax.set_ylim3d(ys.min(),ys.max())
        return ax
    #}}}

#{{{ structured array helper functions
def make_bar_graph_indices(mystructarray,list_of_text_fields,
        recursion_depth = 0,
        verbose = False,
        spacing = 0.1):
    r"This is a recursive function that is used as part of textlabel_bargraph; it does NOT work without the sorting given at the beginning of that function"
    #{{{ if there are still text fields left, then break down the array further, otherwise, just return the indices for this subarray
    if len(list_of_text_fields) > 0:
        unique_values = unique(mystructarray[list_of_text_fields[0]])# the return_index argument doesn't do what it's supposed to all the time, so I have to manually find the start indices, as given in the following line
        start_indices = [nonzero(mystructarray[list_of_text_fields[0]] == val)[0][0] for val in unique_values]
        # find the structured array for the unique value
        index_values = []
        label_values = []
        start_indices = r_[start_indices,len(mystructarray)] # I add this so I can do the next step
        if verbose: print 'recursion depth is',recursion_depth,'and I am analyzing',list_of_text_fields[0],': ',
        if verbose: print 'I found these unique values:',unique_values,'at these start indices:',start_indices[:-1]
        for k in range(0,len(start_indices)-1):
            if verbose: print 'recursion depth is',recursion_depth,'and I am analyzing',list_of_text_fields[0],': ',
            if verbose: print 'trying to extract unique value',unique_values[k],'using the range',start_indices[k],start_indices[k+1]
            if verbose: print 'which has this data'
            indiv_struct_array = mystructarray[start_indices[k]:start_indices[k+1]]
            if verbose: print lsafen(indiv_struct_array)
            these_index_values,these_labels = make_bar_graph_indices(indiv_struct_array,list_of_text_fields[1:],recursion_depth = recursion_depth+1,verbose = verbose)
            index_values.append(these_index_values)
            label_values.append([str(unique_values[k])+','+j for j in these_labels])
        #{{{ scale the result of each call down to the equal size (regardless of number of elements), shift by the position in this array, and return
        if verbose: print 'recursion depth is',recursion_depth,'and I just COMPLETED THE LOOP, which gives a list of index values like this',index_values
        max_indices = max(array(map(len,index_values),dtype='double'))# the maximum width of the array inside
        index_values = map(lambda x: x+(max_indices-len(x))/2.0,index_values)# if the bar is less than max indices, shift it over, so it's still in the center
        if verbose: print 'recursion depth is',recursion_depth,'and I centered each set like this',index_values
        index_values = map(lambda x: x/max_indices*(1-spacing)+(1-spacing)/2,index_values)# scale down, so the width from left edge of bar to right edge of largest bar runs 0--> 1
        if verbose: print 'recursion depth is',recursion_depth,'and I scaled down so each runs zero to one*(1-spacing) (centered) like this',index_values
        # this adds an index value, and also collapses down to a single dimension list
        retval_indices = [x+num for num,val in enumerate(index_values) for x in val]
        # now collapse labels down to a single dimension
        retval_labels = [k for j in label_values for k in j]
        if verbose: print 'recursion depth is',recursion_depth,'and I am passing up indices',retval_indices,'and labels',retval_labels
        return retval_indices,retval_labels
        #}}}
    else:
        if verbose: print 'recursion depth is',recursion_depth,
        N = len(mystructarray)
        if verbose: print 'hit innermost (no text labels left) and passing up a list of indices that looks like this:',r_[0:N]
        return r_[0:N],['']*N
    #}}}

def textlabel_bargraph(mystructarray,othersort = None,spacing = 0.1,verbose = False,ax = None,tickfontsize = 8):
    if ax is None:
        thisfig = gcf()
        ax = thisfig.add_axes([0.2,0.5,0.8,0.5])
        try:
            ax.tick_params(axis = 'both',which = 'major',labelsize = tickfontsize)
            ax.tick_params(axis = 'both',which = 'minor',labelsize = tickfontsize)
        except:
            print 'Warning, in this version I can\'t set the tick params method for the axis'
    #{{{ find the text fields, put them first, and sort by them
    mystructarray = mystructarray.copy()
    list_of_text_fields = [str(j[0]) for j in mystructarray.dtype.descr if j[1][0:2] == '|S']
    mystructarray = mystructarray[list_of_text_fields + [x[0]
        for x in mystructarray.dtype.descr
        if x[0] not in list_of_text_fields]]
    mystructarray.sort()
    if verbose: print 'test --> now, it has this form:',lsafen(mystructarray)
    #}}}
    error_fields = [str(j) for j in mystructarray.dtype.names if j[-6:] == '_ERROR']
    if len(error_fields) > 0:
        mystructarray_errors = mystructarray[error_fields]
        if verbose: "found error fields:",mystructarray_errors
    mystructarray = mystructarray[[str(j) for j in mystructarray.dtype.names if j not in error_fields]]
    if othersort is not None:
        list_of_text_fields.append(othersort)
    if verbose: print 'list of text fields is',lsafen(list_of_text_fields)
    indices,labels = make_bar_graph_indices(mystructarray,list_of_text_fields,verbose = verbose,spacing = spacing)
    temp = zip(indices,labels)
    if verbose: print '(indices,labels) (len %d):'%len(temp),lsafen(temp)
    if verbose: print 'I get these labels (len %d):'%len(labels),labels,'for the data (len %d)'%len(mystructarray),lsafen(mystructarray)
    indices = array(indices)
    indiv_width = min(diff(indices))*(1-spacing)
    remaining_fields = [x for x in mystructarray.dtype.names if x not in list_of_text_fields] # so they are in the right order, since set does not preserve order
    if verbose: print 'The list of remaining (i.e. non-text) fields is',lsafen(remaining_fields)
    colors = ['b','g','r','c','m','k']
    rects = []
    for j,thisfield in enumerate(remaining_fields):
        field_bar_width = indiv_width/len(remaining_fields)
        thiserror = None
        if thisfield+'_ERROR' in error_fields:
            thiserror = mystructarray_errors[thisfield+'_ERROR']
        try:
            rects.append(ax.bar(indices+j*field_bar_width,
                    mystructarray[thisfield],
                    field_bar_width,color = colors[j],
                    yerr = thiserror,#just to test
                    ecolor = 'k',
                    label = '$%s$'%thisfield))
        except:
            raise CustomError('Problem with bar graph: there are %d indices, but %d pieces of data'%(len(indices),len(mystructarray[thisfield])),'indices:',indices,'data',mystructarray[thisfield])
    ax.set_xticks(indices+indiv_width/2)
    ax.set_xticklabels(labels)
    ax.legend([j[0] for j in rects],
            ['$%s$'%j for j in remaining_fields],loc = 'best')
    return

#{{{ a better version?
def othergridandtick(ax,rotation=(0,0),precision=(2,2),labelstring=('',''),gridcolor=r_[0,0,0],y = True,x = True,spines = None):
    #{{{ taken from matplotlib examples
    def adjust_spines(ax,spines):
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position(('outward',5)) # outward by 5 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([],minor = False)
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([],minor = False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #}}}
    if spines is not None:
        adjust_spines(gca(),spines = spines)
    if x:
        #{{{x ticks
        # determine the size
        ax.xaxis.set_major_locator(MaxNLocator(10)) # could use multiplelocator if it keeps try to do multiples of 2
        ax.xaxis.set_minor_locator(MaxNLocator(50))
        #}}}
    if y:
        #{{{ y ticks
        ax.yaxis.set_major_locator(MaxNLocator(10))
        ax.yaxis.set_minor_locator(MaxNLocator(50))
        #}}}
    grid(True,which='major',color=gridcolor,alpha=0.2,linestyle='-')
    grid(True,which='minor',color=gridcolor,alpha=0.1,linestyle='-')
    if x:
        labels = ax.get_xticklabels()
        setp(labels,rotation=rotation[0],fontsize=10)
    if y:
        labels = ax.get_yticklabels()
        setp(labels,rotation=rotation[1],fontsize=10)
    return
#}}}
#{{{ old grid and tick
def gridandtick(ax,rotation=(0,0),precision=(2,2),
        labelstring=('',''),gridcolor=r_[0,0,0],
        formatonly = False,fixed_y_locator = None,
        logarithmic = False,use_grid = True,
        spines = None,y = True):
    #{{{ taken from matplotlib examples
    def adjust_spines(ax,spines):
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position(('outward',5)) # outward by 5 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine
        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([],minor = False)
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([],minor = False)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    #}}}
    if spines is not None:
        adjust_spines(gca(),spines = spines)
    if not formatonly:
        #{{{x ticks
        # determine the size
        width = abs(diff(ax.get_xlim()))
        if width==0:
            raise CustomError('x axis width is zero')
        widthexp = floor(log(width)/log(10.))-1
        scalefactor = 10**widthexp
        width /= scalefactor
        majorLocator = MultipleLocator(5*scalefactor)
        #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[0]+'f'+labelstring[0])# labelstring can be used, for instance, for pi
        #ax.xaxis.set_major_formatter(majorFormatter)
        minorLocator   = MultipleLocator(1*scalefactor)
        ax.xaxis.set_major_locator(majorLocator)
        #for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)
        #}}}
        if y:
            #{{{ y ticks
            width = abs(diff(ax.get_ylim()))
            if width==0:
                raise CustomError('y axis width is zero')
            widthexp = floor(log(width)/log(10.))-1
            scalefactor = 10**widthexp
            width /= scalefactor
            if fixed_y_locator == None:
                if logarithmic:
                    majorLocator = LogLocator(10)
                else:
                    majorLocator   = MultipleLocator(5*scalefactor)
            else:
                majorLocator   = MultipleLocator(fixed_y_locator[4::5])
            #majorFormatter = FormatStrFormatter('%0.'+'%d'%precision[1]+'f'+labelstring[1])# labelstring can be used, for instance, for pi
            #ax.yaxis.set_major_formatter(majorFormatter)
            if fixed_y_locator == None:
                if logarithmic:
                    minorLocator = LogLocator(10,subs=r_[0:11])
                else:
                    minorLocator   = MultipleLocator(1*scalefactor)
            else:
                minorLocator   = FixedLocator(fixed_y_locator)
            ax.yaxis.set_major_locator(majorLocator)
            #for the minor ticks, use no labels; default NullFormatter
            ax.yaxis.set_minor_locator(minorLocator)
            #}}}
    grid(use_grid,which='major',color=gridcolor,alpha=0.15,linestyle='-')
    grid(use_grid,which='minor',color=gridcolor,alpha=0.125,linestyle='-')
    labels = ax.get_xticklabels()
    setp(labels,rotation=rotation[0],fontsize=10)
    if y:
        labels = ax.get_yticklabels()
        setp(labels,rotation=rotation[1],fontsize=10)
    fig = gcf()
    fig.autofmt_xdate()
    return
def gridon(gridcolor=r_[0,0,0]):
    grid(True,which='major',color=gridcolor,alpha=0.1,linestyle='-')
    grid(True,which='minor',color=gridcolor,alpha=0.05,linestyle='-')
#}}}

#{{{ plot wrapper
global OLDplot
OLDplot = plot
global myplotfunc
myplotfunc = OLDplot
def whereblocks(a): # returns contiguous chunks where the condition is true
    parselist = where(a)[0]
    jumps_at = where(diff(parselist)>1)[0]+1
    retlist = []
    lastjump = 0
    for jump in jumps_at:
        retlist += [parselist[lastjump:jump]]
        lastjump = jump
    retlist += [parselist[lastjump:]]
    return retlist
def autolegend(*args,**kwargs):
    #lg = legend(legendstr,'best'),loc = 2, borderaxespad = 0.)
    match_colors = False
    if 'match_colors' in kwargs.keys():
        match_colors = kwargs.pop('match_colors')
    alpha = 0.45
    if 'alpha' in kwargs.keys():
        alpha = kwargs.pop('alpha')
    if 'ax' in kwargs.keys():
        ax_list = [kwargs.pop('ax')]
    else:
        ax_list = [gca()]
    if 'ax2' in kwargs.keys():
        ax_list.append(kwargs.pop('ax2'))
    for ax in ax_list:
        if len(args)==0:
            lg = ax.legend(**kwargs)
        elif len(args)==1:
            lg = ax.legend(args[0],**kwargs)
        else:
            lg = ax.legend(args[0],args[1],**kwargs)
        if lg is None:
            raise ValueError("Warning! you called autolegend, but you don't seem to have anything labeled!!")
        else:
            lg.get_frame().set_alpha(alpha)
    if match_colors:
        for line, txt in zip(lg.get_lines(), lg.get_texts()): # from http://stackoverflow.com/questions/13828246/matplotlib-text-color-code-in-the-legend-instead-of-a-line 
                    txt.set_color(line.get_color())  
                    txt.set_alpha(line.get_alpha())  
    return lg
def autopad_figure(pad = 0.2,centered = False):
    #{{{ solve the axis issue --> this does just the left
    fig = gcf()
    ax = gca()
    labelsets = [] 
    #labelsets.append(('left',ax.get_yticklabels()))
    #labelsets.append(('left',ax.get_yticklines()))
    #labelsets.append(('right',ax.get_yticklines()))
    labelsets.append(('left',[ylabel(ax.get_ylabel())]))
    #labelsets.append(('bottom',ax.get_xticklabels()))
    #labelsets.append(('bottom',ax.get_xticklines()))
    if len(ax.get_xlabel()) > 0:
        labelsets.append(('bottom',[xlabel(ax.get_xlabel())]))
    #labelsets.append(('top',ax.get_xticklines()))
    if len(ax.get_title()) > 0:
        pass #labelsets.append(('top',[title(ax.get_title())]))
    compto = {}
    def on_draw(event):
       # find the sum of the widths of all things labeled with a 'y'
       spkwargs = {}
       compto['bottom'] = fig.subplotpars.bottom
       compto['left'] = fig.subplotpars.left
       compto['right'] = fig.subplotpars.right
       compto['top'] = fig.subplotpars.top
       for axisn in ['left','bottom','top','right']:
           bboxes = []
           labellist = [x[1] for x in labelsets if x[0] is axisn]
           for labels in labellist:
               for label in labels:
                   if type(label) is Line2D:
                       pass # just rely on the pad
                       #if any(map(lambda x: x == label.get_transform(),[ax.transData,ax.transAxes,fig.transFigure,None])):
                       #    print 'found it'
                       #else:
                       #    print 'didn not find it'
                       #bbox = label.get_window_extent(fig.canvas).inverse_transformed(ax.transData).inverse_transformed(fig.transFigure)
                   else:
                       try:
                           bbox = label.get_window_extent()
                       except:
                           raise CustomError('type of label = ',type(label))
                   # the figure transform goes from relative coords->pixels and we
                   # want the inverse of that
                   bboxes.append(bbox)
               # this is the bbox that bounds all the bboxes, again in relative
               # figure coords
           l = 0 
           if len(labellist):
               bbox = mtransforms.Bbox.union(bboxes)
               bboxi = bbox.inverse_transformed(fig.transFigure)
               if axisn in ['left','right']:
                   l = bboxi.width
               if axisn in ['top','bottom']:
                   l = bboxi.height
           l += pad
           if axisn in ['top','right']:
               l = 1-l
               if compto[axisn] > l:
                   spkwargs.update({axisn:l})
           else:
               if compto[axisn] < l:
                   spkwargs.update({axisn:l})
       try:
           if len(spkwargs) > 0:
               if centered and 'left' in spkwargs.keys() and 'right' in spkwargs.keys():
                   big = max(r_[spkwargs['left'],1-spkwargs['right']])
                   spkwargs.update({'left':big,'right':1-big})
               fig.subplots_adjust(**spkwargs) # pad a little
               #print "adjusted to",spkwargs
               fig.canvas.draw()
       except:
           raise CustomError('spwargs = ',spkwargs)
       return False
    fig.canvas.mpl_connect('draw_event', on_draw)
    fig.subplots_adjust(left = 0, right = 1, top = 1, bottom =0)
    fig.canvas.draw()
    #}}}
def expand_x(*args):
    # this is matplotlib code to expand the x axis
    ax = gca()
    xlims = array(ax.get_xlim())
    xlims.sort()
    width = abs(diff(xlims))
    xlims[0] -= width/10
    xlims[1] += width/10
    if len(args) > 0:
        if len(args) == 1 and type(args) is tuple:
            args = args[0]
        xlims = [xlims[j] if args[j] is None else args[j] for j in range(0,2)]
    ax.set_xlim(xlims)
def expand_y(*args):
    # this is matplotlib code to expand the x axis
    ax = gca()
    ylims = array(ax.get_ylim())
    width = abs(diff(ylims))
    ylims[0] -= width/10
    ylims[1] += width/10
    if len(args) > 0:
        if len(args) == 1 and type(args) is tuple:
            args = args[0]
        ylims = [ylims[j] if args[j] is None else args[j] for j in range(0,2)]
    ax.set_ylim(ylims)
def plot_label_points(x,y,labels,**kwargs_passed):
    kwargs = {'alpha':0.5,'color':'g','ha':'left','va':'center','rotation':0,'size':14}
    kwargs.update(kwargs_passed)
    for j in range(0,len(labels)):
        text(x[j],y[j],labels[j],**kwargs)
def addlabels(labelstring,x,y,labels):
    r'obsolete -- use plot_label_points'
    for j in range(0,len(labels)):
        text(x[j],y[j],labelstring%labels[j],alpha=0.5,color='g',ha='left',va='top',rotation=0)
def plot_color_counter(*args,**kwargs):
    if 'ax' in kwargs.keys():
        ax = kwargs.pop('ax')
    else:
        ax = gca()
    if len(args)>0:
        try:
            ax._get_lines.count = args[0] # set the value of the color counter
        except:
            ax._get_lines.color_cycle = args[0] # set the value of the color counter
    try: # this is different depending on the version of matlablike
        retval = ax._get_lines.count
    except:
        retval = ax._get_lines.color_cycle
    return retval
def contour_plot(xvals,yvals,zvals,color = 'k',alpha = 1.0,npts = 300,**kwargs):
    if 'inline_spacing' in kwargs.keys():
        inline_spacing = kwargs.pop('inline_spacing')
    else:
        inline_spacing = 20
    xi = linspace(xvals.min(),xvals.max(),npts)
    yi = linspace(yvals.min(),yvals.max(),npts)
    #{{{ show the diffusivity
    #plot(array(xvals),array(yvals),'k')# to show where everything is
    zi = scipy_griddata((xvals,yvals),
        zvals,
        (xi[None,:],yi[:,None]))
    zi_min = zi[isfinite(zi)].min()
    zi_max = zi[isfinite(zi)].max()
    levels = r_[zi_min:zi_max:40j]
    CS = contour(xi,yi,zi,levels,colors = color,
            alpha = 0.25*alpha)
    oldspacing = levels[1]-levels[0]
    levels = r_[zi_min:zi_max:oldspacing*5]
    try:
        CS = contour(xi,yi,zi,levels,colors = color,
            alpha = alpha,**kwargs)
    except:
        raise CustomError("Is there something wrong with your levels?:",levels,"min z",zi_min,"max z",zi_max)
    clabel(CS,fontsize = 9,inline = 1,
        #fmt = r'$k_\sigma/k_{\sigma,bulk} = %0.2f$',
        fmt = r'%0.2f',
        use_clabeltext = True,
        inline_spacing = inline_spacing,
        alpha = alpha)
    #}}}
def giveSpace(spaceVal = 0.1):#{{{
    """
    This should be called after plot(data) call. This function will put a white spacing around any data and associated error on both x and y dimensions
    """
    ax = gca()
    xticks = ax.get_xticks()
    xspace = xticks.max()*spaceVal
    xmin = xticks.min() - xspace
    xmax = xticks.max() + xspace
    ax.set_xlim(xmin,xmax)
    # y axis 
    yticks = ax.get_yticks()
    yspace = yticks.max()*spaceVal
    ymin = yticks.min() - yspace
    ymax = yticks.max() + yspace
    ax.set_ylim(ymin,ymax)#}}}

def plot_updown(data,axis,color1,color2,symbol = '',**kwargs):
    if symbol == '':
        symbol = 'o'
    change = r_[1,diff(data.getaxis(axis))]
    changemask = change > 0
    if 'force_color' in kwargs.keys() and kwargs['force_color'] == True:
        if hasattr(data,'other_info'):
            if 'plot_color' in data.other_info.keys():
                data.other_info.pop('plot_color')
    plot(data[axis,changemask],color1+symbol,**kwargs)
    if len(kwargs) > 0 and 'label' in kwargs.keys(): kwargs.pop('label') # if I'm doing a legend, I want it on the first
    plot(data[axis,~changemask],color2+symbol,**kwargs)
    return
def nextfigure(figurelist,name):
    'obsolete -- now use class'
    if isinstance(figurelist,figlist_var):
        figurelist.next(name)
        return figurelist
    else:
        print 'Boo! not a new style name!'
    verbose = False # a good way to debug
    if verbose: print lsafe('DEBUG figurelist, called with',name)
    if name in figurelist:
        fig = figure(figurelist.index(name)+1)
        if verbose: print lsafen('in',figurelist,'at figure',figurelist.index(name)+1,'switched figures')
    else:
        fig = figure(len(figurelist)+1)
        fig.add_subplot(111)
        if verbose: print lsafen('added, figure',len(figurelist)+1,'because not in figurelist',figurelist)
        figurelist.append(name)
    return figurelist
def figlistret(first_figure,figurelist,*args,**kwargs):
    if first_figure is None:
        raise ValueError("Don't know how to do this when you run interactively")
    else:
        args += (figurelist,)
        if len(args) == 1:
            return args[0]
        else:
            return args
def figlistini_old(first_figure):
    if isinstance(first_figure,figlist_var):
        return first_figure
    verbose = False
    if verbose: print lsafe('DEBUG: initialize figlist')
    if first_figure == None:
        if verbose: print lsafen('empty')
        return []
    else:
        if verbose: print lsafen(first_figure.figurelist)
        return first_figure
class figlist(object):
    def __init__(self,*arg,**kwargs):
        self.verbose = False
        self.black = 0.9
        self.env = ''
        if 'mlab' in kwargs.keys():
            self.mlab = kwargs.pop('mlab')
        if 'env' in kwargs.keys():
            self.env = kwargs.pop('env')
        if 'verbose' in kwargs.keys():
            self.verbose = kwargs.pop('verbose')
        if self.verbose: print lsafe('DEBUG: initialize figlist')
        if len(arg) == 0:
            if self.verbose: print lsafen('empty')
            self.figurelist = []
        else:
            if self.verbose: print lsafen(arg[0])
            self.figurelist = arg[0]
        if len(kwargs) > 0:
            self.figurelist.append(kwargs)
        self.units = {}
        self.autolegend_list = {}
        self.twinx_list = {}
        self.basename = None
        return
    def twinx(self,autopad = True,orig = False,color = None):
        #self.figurelist.insert(self.figurelist.index(self.current),{'autopad':False}) #doesn't work because it changes the figure number; I can get the number with fig = gcf(); fig.number, but I can't set it; it would be best to switch to using a list that contains all the figure numbers to match all their names -- or alternatively, note that matplotlib allows you to give them names, though I don't know how that works
        if self.current in self.twinx_list.keys():
            ax1,ax2 = self.twinx_list[self.current]
        else:
            autopad_figure()
            ax1 = gca()
            twinx()
            ax2 = gca()
            self.twinx_list[self.current] = (ax1,ax2)
            if color is not None:
                print ax2.spines
                ax2.tick_params(axis = 'y',colors = color)
                ax2.yaxis.label.set_color(color)
                ax2.spines['right'].set_color(color)
        if orig:
            sca(ax1)
            return ax1
        else:
            sca(ax2)
            return ax2
    def use_autolegend(self,value = None):
        'No argument sets to true if it\'s not already set'
        if value is None:
            if not self.current in self.autolegend_list.keys():
                self.autolegend_list.update({self.current:True})
            else: #leave it alone
                return
        else: #passed an explicit value
            self.autolegend_list.update({self.current:value})
            return
    def push_marker(self):
        if not hasattr(self,'pushlist'):
            self.pushlist = []
        self.pushlist.append(self.current)
        return
    def pop_marker(self):
        self.next(self.pushlist.pop())
        return
    def next(self,name,legend = False,boundaries = None,twinx = None,**kwargs):
        if self.basename is not None: #basename for groups of figures
            name = self.basename + '_' + name
        if name.find('/') > 0:
            raise ValueError("don't include slashes in the figure name, that's just too confusing")
        if self.verbose: print lsafe('DEBUG figurelist, called with',name)
        if name in self.figurelist:
            if hasattr(self,'mlab'):
                fig = mlab.figure(self.figurelist.index(name)+1,bgcolor = (1,1,1),**kwargs)
                fig.scene.render_window.aa_frames = 20
                fig.scene.anti_aliasing_frames = 20
            else:
                fig = figure(self.figurelist.index(name)+1,**kwargs)
            self.current = name
            if self.verbose: print lsafen('in',self.figurelist,'at figure',self.figurelist.index(name)+1,'switched figures')
        else:
            self.current = name
            if boundaries == False:
                self.setprops(boundaries = False)
            if legend:
                if 'figsize' not in kwargs.keys():
                    kwargs.update({'figsize':(12,6)})
                if hasattr(self,'mlab'):
                    fig = mlab.figure(len(self.figurelist)+1,bgcolor = (1,1,1),**kwargs)
                    fig.scene.render_window.aa_frames = 8
                else:
                    fig = figure(len(self.figurelist)+1,**kwargs)
                fig.add_axes([0.075,0.2,0.7,0.7]) # l b w h
                self.use_autolegend('outside')
            else:
                fig = figure(len(self.figurelist)+1,**kwargs)
                if hasattr(self,'mlab'):
                    fig = self.mlab.figure(len(self.figurelist)+1,bgcolor = (1,1,1),**kwargs)
                    fig.scene.render_window.aa_frames = 8
                else:
                    fig = figure(len(self.figurelist)+1,**kwargs)
                if twinx is not None:
                    fig.add_subplot(111)
            if self.verbose: print lsafen('added, figure',len(self.figurelist)+1,'because not in figurelist',self.figurelist)
            self.figurelist.append(name)
        if twinx is not None:
            if twinx == 0:
                self.twinx(orig = True)
                fig = gcf()
            elif twinx == 1:
                self.twinx()
                fig = gcf()
            else:
                raise ValueError('If you pass twinx, pass 0 for the original or 1 for the right side')
        return fig
    def plot(self,*args,**kwargs):
        if 'label' in kwargs.keys():
            self.use_autolegend()
        human_units = True
        if 'human_units' in kwargs.keys():
            human_units = kwargs.pop('human_units')
        if human_units:
            firstarg = self.check_units(args[0],0,1) # check units, and if need be convert to human units, where x is the first dimension and y is the last
        else:
            firstarg = args[0]
        if 'label' not in kwargs.keys() and isinstance(args[0],nddata):
            thisname = args[0].name()
            if thisname is not None:
                kwargs['label'] = thisname
        retval = plot(*tuple((firstarg,)+args[1:]),**kwargs)#just a placeholder for now, will later keep units + such
        ax = gca()
        if ax.get_title() is None or len(ax.get_title()) == 0:
            try:
                title(self.current)
            except:
                title('untitled')
        return retval
    def check_units(self,testdata,x_index,y_index):
        if isinstance(testdata,nddata):
            testdata = testdata.copy().human_units()
            if len(testdata.dimlabels) > 1:
                if not hasattr(self,'current'):
                    raise ValueError("give your plot a name (using .next()) first! (this is used for naming the PDF's etc)")
                if self.current in self.units.keys():
                        theseunits = (testdata.get_units(testdata.dimlabels[x_index]),testdata.get_units(testdata.dimlabels[y_index]))
                        if theseunits != self.units[self.current] and theseunits[0] != self.units[self.current]:
                                raise ValueError("the units don't match (old units %s and new units %s)! Figure out a way to deal with this!"%(theseunits,self.units[self.current]))
                else:
                    if isinstance(testdata,nddata):
                        self.units[self.current] = (testdata.get_units(testdata.dimlabels[x_index]),testdata.get_units(testdata.dimlabels[y_index]))
            else:
                if not hasattr(self,'current'):
                    self.next('default')
                if self.current in self.units.keys():
                    theseunits = (testdata.get_units(testdata.dimlabels[x_index]))
                    testunits = self.units[self.current]
                    if theseunits != testunits:
                        if type(testunits) is tuple and testunits[1] is None:
                            pass
                        else:
                            raise ValueError("the units don't match (old units %s and new units %s)! Figure out a way to deal with this!"%(theseunits,self.units[self.current]))
                else:
                    self.units[self.current] = (testdata.get_units(testdata.dimlabels[x_index]))
        return testdata
    def adjust_spines(self,spines):
        ax = gca()
        #{{{ taken from matplotlib examples
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position(('outward',10)) # outward by 10 points
                spine.set_smart_bounds(True)
            else:
                spine.set_color('none') # don't draw spine

        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([])

        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([])
        #}}}
    def grid(self):
        ax = gca()
        if self.black:
            gridandtick(ax,gridcolor = r_[0.5,0.5,0.5])
        else:
            gridandtick(ax,gridcolor = r_[0,0,0])
        return
    def image(self,*args,**kwargs):
        firstarg = self.check_units(args[0],-1,0) # check units, and if need be convert to human units, where x is the last dimension and y is the first
        if self.black and 'black' not in kwargs.keys():
            kwargs.update({'black':self.black})
        retval = image(*tuple((firstarg,)+args[1:]),**kwargs)#just a placeholder for now, will later keep units + such
        ax = gca()
        if ax.get_title() is None or len(ax.get_title()) == 0:
            title(self.current)
        return retval
    def text(self,mytext):
        self.setprops(print_string = mytext)
    def setprops(self,**kwargs):
        self.figurelist.append(kwargs)
    def show_prep(self):
        for k,v in self.autolegend_list.items():
            kwargs = {}
            if v:
                if type(v) is str:
                    if v[0:7] == 'colored':
                        kwargs.update(dict(match_colors = True))
                        v = v[7:]
                        if v == '':
                            v = True
                    if v == 'outside':
                        kwargs.update(dict(bbox_to_anchor=(1.05,1),loc = 2,borderaxespad=0.))
                self.next(k)
                try:
                    autolegend(**kwargs)
                except:
                    try:
                        self.twinx(orig = True)
                        autolegend(**kwargs)
                    except:
                        raise CustomError('error while trying to run autolegend for',k)
    def show(self,*args,**kwargs):
        verbose = False
        if 'verbose' in kwargs.keys():
            verbose = kwargs.pop('verbose')
        if len(kwargs) > 0:
            raise ValueError("didn't understand kwargs "+repr(kwargs))
        self.show_prep()
        #{{{ just copy from fornnotebook to get the print string functionality
        kwargs = {}
        for j,figname in enumerate(self.figurelist):
            if verbose: print "showing figure"+lsafen(figname)
            if type(figname) is dict:
                kwargs.update(figname)
                if 'print_string' in kwargs:
                    print '\n\n'
                    print kwargs.pop('print_string')
                    print '\n\n'
        #}}}
        if len(args) == 1:
            if (args[0][:-4] == '.pdf') or (args[0][:-4] == '.png') or (args[0][:-4] == '.jpg'):
                print "you passed me a filename, but I'm just burning it"
        if hasattr(self,'mlab'):
            print "running mlab show!"
            self.mlab.show()
        else:
            #print "not running mlab show!"
            show()
    def header(self,number_above,input_string):
        header_list = ['\\section','\\subsection','\\subsubsection','\\paragraph','\\subparagraph']
        self.text(header_list[number_above+1]+'{%s}'%input_string)
        return number_above + 1
def text_on_plot(x,y,thistext,coord = 'axes',**kwargs):
    ax = gca()
    if coord == 'axes':
        newkwargs = {'transform':ax.transAxes,'size':'x-large',"horizontalalignment":'center'}
    elif coord == 'data':
        print "Yes, I am using data transform"
        newkwargs = {'transform':ax.transData,'size':'small',"horizontalalignment":'right'}
    color = None
    if 'match_data' in kwargs.keys():
        if type(kwargs['match_data']) is list:
            color = kwargs['match_data'][-1].get_color() # get the color of the last line
        elif kwargs['match_data'].get_plot_color() is not None:
            color = kwargs['match_data'].get_plot_color() # don't know when this works, but apparently, it does!
        if color is not None:
            newkwargs.update({'color':color})
        else:
            raise CustomError('You passed match_data to text_on_plot, but I can\'t find a color in the object')
        kwargs.pop('match_data')
    newkwargs.update(kwargs)
    return text(x,y,thistext,**newkwargs)
def unitify_axis(myy,axis_name,is_axis = True):
    'this just generates an axis label with appropriate units'
    if is_axis:
        yunits = myy.units_texsafe(axis_name)
        j = axis_name.find('_')
        if j > -1:
            prevword = axis_name[0:j]
            if j+1< len(axis_name):
                followword = axis_name[j+1:]
            else:
                followword = []
            k = followword.find(' ')
            if k > -1 and k < len(followword):
                followword = followword[:k]
            k = followword.find('_')
            if len(followword) > 0:
                if not (k > -1) and (len(prevword) < 2 or len(followword) < 2):
                    if len(followword) > 1:
                        axis_name = axis_name[:j+1+len(followword)]  + '}$' + axis_name[j+1+len(followword):]
                        axis_name = axis_name[:j+1] + '{' + axis_name[j+1:]
                    else:
                        axis_name = axis_name[0:j+2] + '$' + axis_name[j+2:]
                    axis_name = '$'+axis_name
        if myy.get_prop('FT'):
            axis_name = r'F{'+axis_name+r'}'
    else:
        yunits = myy.units_texsafe()
    if yunits is not None:
        axis_name = axis_name + ' / ' + yunits
    return axis_name
def plot(*args,**kwargs):
    global myplotfunc
    has_labels = False
    #{{{ deal with axes and some other kwargs
    if 'ax' in kwargs.keys():
        ax = kwargs.pop('ax')
    else:
        ax = gca()
    human_units = False
    if 'human_units' in kwargs.keys():
        human_units = kwargs.pop('human_units')
    label_format_string = None
    if 'label_format_string' in kwargs.keys():
        label_format_string = kwargs.pop('label_format_string')
    normalize = False
    if 'normalize' in kwargs.keys():
        normalize = kwargs.pop('normalize')
    #}}}
    myplotfunc = ax.plot # default
    #{{{ all possible properties
    myformat = None 
    myxlabel = None
    myylabel = None
    myx = None
    myy = None
    #}}}
    #{{{assign all the possible combinations
    if len(args)==1:
        myy = args[0]
    elif (len(args)==2) and (type(args[1]) is str):
        myy = args[0]
        myformat = args[1]
    else:
        myx = args[0]
        myy = args[1]
    if len(args)==3:
        myformat = args[2]
    if isscalar(myx):
        myx = array([myx])
    if isscalar(myy):
        myy = array([myy])
    #}}}
    #{{{ parse nddata
    if isinstance(myy,nddata):
        myy = myy.copy()
        if human_units: myy = myy.human_units()
        if myy.get_plot_color() is not None\
            and 'color' not in kwargs.keys():# allow override
            kwargs.update({'color':myy.get_plot_color()})
        if myy.name() is not None:
            myylabel = myy.name()
        else:
            myylabel = 'data'
        myylabel = unitify_axis(myy,myylabel,is_axis = False)
        if (len(myy.dimlabels)>0):
            myxlabel = myy.dimlabels[0]
            xunits = myy.units_texsafe(myxlabel)
            if xunits is not None:
                myxlabel += ' / ' + xunits
        if (myx == None):
            try:
                myx = myy.getaxis(myy.dimlabels[0])
            except:
                myx = r_[0:myy.data.shape[0]]
        if type(myy.data_error) is ndarray and len(myy.data_error)>0: #then this should be an errorbar plot
            def thiserrbarplot(*tebargs,**tebkwargs):
                if type(tebargs[-1]) is str:
                    tebkwargs.update({'fmt':tebargs[-1]})
                    ax.errorbar(*tebargs[:-1],**tebkwargs)
                else:
                    ax.errorbar(*tebargs,**tebkwargs)
            myplotfunc = thiserrbarplot
            #{{{ pop any singleton dims
            myyerror = myy.get_error()
            myyerror = squeeze(myyerror)
            #}}}
            kwargs.update({'yerr':myyerror})
            valueforxerr = myy.get_error(myy.dimlabels[0])
            if valueforxerr != None: # if we have x errorbars too
                #print "DEBUG decided to assign to xerr:",valueforxerr
                kwargs.update({'xerr':valueforxerr})
        #{{{ deal with axis labels along y
        try:
            yaxislabels = myy.getaxis(myy.dimlabels[-1])
        except:
            pass
        # at this point, if there is no axis label, it will break and go to pass
        if yaxislabels is not None:
            if len(yaxislabels) > 0:
                if type(yaxislabels[0]) is string_:
                    has_labels = True
                elif label_format_string is not None:
                    yaxislabels = [label_format_string%j for j in yaxislabels]
                    has_labels = True
        #}}}
        myy = squeeze(myy.data)
    #}}}
    #{{{ semilog where appropriate
    # I don't want to semilog and I'm not sure why this got changed...
    #if (myx != None) and (len(myx)>1): # by doing this and making myplotfunc global, we preserve the plot style if we want to tack on one point
    #    try:
    #        b = diff(log10(myx))
    #    except:
    #        raise CustomError('likely a problem with the type of the x label, which is',myx)
    #    if (size(b)>3) and all(abs((b-b[0])/b[0])<1e-4) and not ('nosemilog' in kwargs.keys()):
    #        myplotfunc = ax.semilogx
    if ('nosemilog' in kwargs.keys()):
        #print 'this should pop nosemilog'
        kwargs.pop('nosemilog')
    if 'plottype' in kwargs.keys():
        if kwargs['plottype'] == 'semilogy':
            myplotfunc = ax.semilogy
        elif kwargs['plottype'] == 'semilogx':
            myplotfunc = ax.semilogx
        elif kwargs['plottype'] == 'loglog':
            myplotfunc = ax.loglog
        kwargs.pop('plottype')
    #}}}
    #{{{ take care of manual colors
    if myformat != None:
        colorpos = myformat.find('#')
        if  colorpos > -1:
            kwargs.update({'color':myformat[colorpos:colorpos+7]})
            myformat = myformat[0:colorpos] + myformat[colorpos+7:]
        ##kwargs.update({'fmt':myformat})
        linematched = False
        for linestyle in ['-','--','-.',':','None','  ']:
            if myformat.find(linestyle) > -1:
                linematched = True
                myformat.replace(linestyle,'')
                kwargs.update({'linestyle':linestyle})
        for markerlabel in ['o','.','d']:
            if myformat.find(markerlabel) > -1:
                if not linematched: kwargs.update({'linestyle':''})
                myformat.replace(markerlabel,'')
                kwargs.update({'marker':markerlabel})
        if len(myformat) == 0:
            myformat = None
    #}}}
    if normalize:
        myy /= myy.max()
    #{{{ hsv plots when we have multiple lines
    if len(shape(myy))>1 and sum(array(shape(myy))>1):
        #{{{ hsv plots
        hold(True)
        retval = []
        for j in range(0,myy.shape[1]):
            #{{{ this is the way to assign plot arguments
            plotargs = [myx,myy[:,j],myformat]
            while None in plotargs:
                plotargs.remove(None)
            #}}}
            #{{{ here, i update the kwargs to include the specific color for this line
            newkwargs = kwargs.copy() # kwargs is a dict
            newkwargs.update({'color':cm.hsv(double(j)/double(myy.shape[1]))})
            #}}}
            #{{{ here, I update to use the labels
            if has_labels:
                newkwargs.update({'label':yaxislabels[j]})
            #}}}
            if any(isinf(myy)):
                myy[isinf(myy)] = NaN # added this to prevent an overflow error
            try:
                retval += [myplotfunc(*tuple(plotargs),**newkwargs)]
                #print "\n\n\\begin{verbatim}DEBUG plot:",plotargs,'\nkwargs:\n',newkwargs,'\\end{verbatim}'
            except: 
                raise CustomError("Error trying to plot using function",myplotfunc,len(plotargs),"arguments",plotargs,"of len",map(len,plotargs),"and",len(newkwargs),"options",newkwargs,"of len",map(len,newkwargs.values()))
        #hold(False)
        #}}}
        #}}}
    else:
        plotargs = [myx,real(myy),myformat]
        while None in plotargs:
            plotargs.remove(None)
        try:
            #print 'DEBUG plotting with args',plotargs,'and kwargs',kwargs,'\n\n'
            retval = myplotfunc(*plotargs,**kwargs)
        except:
            raise CustomError('error trying to plot',myplotfunc,
                    '\nlength of the ndarray arguments:',['shape:'+repr(shape(j)) if type(j) is ndarray else j for j in plotargs],
                    '\nsizes of ndarray kwargs',dict([(j,shape(kwargs[j])) if type(kwargs[j]) is ndarray else (j,kwargs[j]) for j in kwargs.keys()]),
                    '\narguments = ',plotargs,
                    '\nkwargs =',kwargs)
    #{{{ attach labels and such
    if (myxlabel!=None):
        ax.set_xlabel(myxlabel)
    if (myylabel!=None):
        ax.set_ylabel(myylabel)
    try:
        ax.axis('tight')
    except:
        raise CustomError('error trying to set axis tight after plot',myplotfunc,'with arguments',plotargs,'and kwargs',kwargs,'\nsizes of arguments:',[shape(j) for j in plotargs],'\nsizes of ndarray kwargs:',dict([(j,shape(kwargs[j])) for j in kwargs.keys() if type(kwargs[j]) is ndarray]))
    #grid(True)
    #}}}
    return retval
#}}}

#{{{subplot_dim
class subplot_dim():
    def __init__(self,firstdim,seconddim):
        self.num = r_[firstdim,seconddim,0]
    def set(self,args,x='',g=True,y='',t='',a=''):
        if type(args) is int:
            number = args
            ax = subplot(*tuple(self.num+r_[0,0,number]))
            xlabel(x)
            ylabel(y)
            title(t)
            grid(g)
        elif (type(args) is tuple) and (len(args) is 3):
            # the second value passed is 
            whichsmall = args[2]
            break_into = args[1]
            number = args[0]
            mydims = self.num*r_[1,break_into,1]+r_[
                    0,0,break_into*(number-1)+whichsmall]
            try:
                ax = subplot(*tuple(mydims))
            except:
                print 'failed trying subplots: ', mydims
                raise
            xlabel(x)
            ylabel(y)
            title(t)
            grid(g)
        else:
            print "problem, need to pass either 1 or 3 arguments to set"
            print 'type of args: ',type(args)
        return ax
#}}}

def spectrogram(waveform,f_start,f_stop,npoints_fdom=40,tdom_div=2):
    #npoints_tdom = int(round(double(waveform.len)/double(npoints_fdom)))*npoints_tdom_mult
    npoints_tdom = waveform.len/tdom_div # this seems to be more legible than above 
    resolution = diff(waveform.x[0:2])

    sigma = abs(f_start-f_stop)/double(npoints_fdom)
    #print "sigma = %f resolution = %f"%(sigma,resolution)
    if sigma<4*resolution:
        sigma = 4*resolution

    waveform.def_filter(sigma,npoints_tdom)# define the filter and number of points for the spectrogram windowing (define the filter such that the points are spaced sigma apart)

    # go through and apply the filter for some range of points

    f_axis = linspace(f_start,f_stop,npoints_fdom)

    specgram = zeros((npoints_fdom,npoints_tdom),dtype="complex128")

    for j in range(0,npoints_fdom):

        t_axis, specgram[j,:] = waveform.do_filter(f_axis[j])
        #plot(t_axis,abs(specgram[j,:])) # leave this in for testing what it does in the fdom
    #image(specgram,y=f_axis/1e6,x=t_axis*1e6) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
    imshow(abs(specgram),extent=(t_axis[0]*1e6,t_axis[-1]*1e6,f_axis[-1]/1e6,f_axis[0]/1e6)) # now do an imagehsv (see if we can make imagerybw) plot of the resulting spectrogram
    return gca()

def image(A,x=[],y=[],**kwargs):
    #{{{ pull out kwargs for imagehsv
    imagehsvkwargs = {}
    spacing = 1
    for k,v in kwargs.items():
        if k in ['black','logscale']:
            imagehsvkwargs[k] = kwargs.pop(k)
        if k in ['spacing']:
            spacing = kwargs.pop(k)
    #}}}
    setlabels = False
    if isinstance(A,nddata):
        setlabels = True
        templabels = list(A.dimlabels)
        x_label = templabels[-1]
        try:
            x = list(A.getaxis(x_label))
        except:
            x = r_[0,ndshape(A)[x_label]]
        x_label = unitify_axis(A,x_label)
        templabels.pop(-1)
        y_label = ''
        if len(templabels) == 1:
            y_label = templabels[0]
            try:
                y = list(A.getaxis(y_label))
            except:
                y = r_[0:A.data.shape[A.axn(y_label)]]
            y_label = unitify_axis(A,y_label)
        else:
            while len(templabels)>0:
                y_label += templabels.pop(0)
                if len(templabels)>0:
                    y_label += '$\\otimes$'
        A = A.data
    if type(x) is list:
        x = array(x)
    if type(y) is list:
        y = array(y)
    if len(x)==0:
        x = [1,A.shape[1]]
    else:
        x = x.flatten()
    if len(y)==0:
        y = [1,A.shape[0]]
    else:
        y = y.flatten()
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    myext = (x[0]-dx/2.,x[-1]+dx/2.,y[-1]+dy/2.,y[0]-dy/2.)
    linecounter = 0
    origAndim = A.ndim
    if A.ndim > 2:
        ax = gca()
        setp(ax.get_yticklabels(),visible = False)
        ax.yaxis.set_ticks_position("none")
    while A.ndim>2:# to substitude for imagehsvm, etc., so that we just need a ersion of ft
        # order according to how it's ordered in the memory
        # the innermost two will form the image -- first add a line to the end of the images we're going to join up
        tempsize = array(A.shape) # make a tuple the right shape
        if linecounter == 0 and spacing < 1.0:
            spacing = round(prod(tempsize[0:-1])) # find the length of the thing not counting the columns
        tempsize[-2] = 2*linecounter + spacing # all dims are the same except the image row, to which I add an increasing number of rows
        #print "iterate (A.ndim=%d) -- now linecounter is "%A.ndim,linecounter
        linecounter += tempsize[-2] # keep track of the extra lines at the end
        A = concatenate((A,nan*zeros(tempsize)),axis=(A.ndim-2)) # concatenate along the rows
        tempsize = r_[A.shape[0:-3],A.shape[-2:]]
        tempsize[-2] *= A.shape[-3]
        A = A.reshape(tempsize) # now join them up
    A = A[:A.shape[0]-linecounter,:] # really I should an extra counter besides linecounter now that I am using "spacing", but leave alone for now, to be sure I don't cute off data
    #line_mask = isnan(A)
    #A[line_mask] = A[logical_not(line_mask)].max()
    #A[line_mask] = 0
    if iscomplex(A).any():
        A = imagehsv(A,**imagehsvkwargs)
        retval = imshow(A,extent=myext,**kwargs)
    else:
        retval = imshow(A,extent=myext,**kwargs)
        colorbar()
    if setlabels:
        xlabel(x_label)
        #print y_label
        ylabel(y_label)
    return retval

def colormap(points,colors,n=256):
    r = interp(linspace(0,1,n),points,colors[:,0].flatten())
    g = interp(linspace(0,1,n),points,colors[:,1].flatten())
    b = interp(linspace(0,1,n),points,colors[:,2].flatten())
    return reshape(r_[r,g,b],(3,n)).T

def imagehsv(A,logscale = False,black = False):
    n = 256
    mask = isnan(A)
    A[mask] = 0
    theta = (n-1.)*mod(angle(A)/pi/2.0,1)# angle in 255*cycles
    hsv = colormap(r_[0.,1./3.,2./3.,1.],double(array([
        [1,0,0],
        [0,1,0],
        [0,0,1],
        [1,0,0]])),n=n)
    hsv_norm = sqrt(sum(hsv*hsv,axis=1))
    hsv_norm = reshape(hsv_norm,(hsv_norm.size,1))
    hsv = hsv/hsv_norm
    colors = hsv[ix_(int32(theta.flatten().round()),[0,1,2])]
    colors = reshape(colors,(A.shape[0],A.shape[1],3))
    mask = mask.reshape(A.shape[0],A.shape[1],1) # reshape the mask into the 3 color shape as well
    mask = tile(mask,(1,3)).reshape(mask.shape[0],mask.shape[1],3) # and copy the mask across all colors
    intensity = abs(A).reshape(A.shape[0],A.shape[1],1)
    intensity /= abs(A).max()
    if logscale:
        intensity = log10(intensity)
    if black:
        if black is True:
            colors = intensity*colors
        else:
            colors = intensity*colors*black + (1.0-black)
        colors[mask] = 1.0
    else:
        colors = 1.0-intensity*(1.0-colors)
        colors[mask] = 0.0
    return colors
