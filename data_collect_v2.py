## Author: Frank Backs
## Email: frank.backs@student.uva.nl
## University of Amsterdam

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy import coordinates as coordfix
from astropy import units as u


small_catalog = []
global small_catalog

###Reading and writing files ---------------------------------------------------

def read(some_filename):
    """
    Reads some file and puts all the lines in a list
    """
    os.chdir("D:/Project data")
    some_file = open(some_filename)
    some_list = []
    for line in some_file:
        some_list.append(line[:-1])
        
    return some_list

def make_file(some_list, filename):
    """
    Makes a file of a list, every entery in the list gets a new line
    """
    os.chdir("D:/Project data")
    some_file = open(filename, 'w')
    for item in some_list:
        some_file.write(str(item) + '\n')
        
    some_file.close()
    
### Catalogs and databases -----------------------------------------------------

def read_database(database):
    """
    Looks for an item at the given coordinates, returns all info about it.
    """
    filehandle = open(database)
    items = [[]]
    for line in filehandle:
        if line[0] != '#':
            if line[0] == '@':
                items.append([])
            else:
                items[-1].append([])
                index = 0
                first = True
                for char in line:
                    if (char != ' ' and char != '\t' and char != '\n') and first: 
                        first_index = index
                        first = False
                    elif (char == ' ' or char == '\t' or char == '\n') and not first:
                        items[-1][-1].append(line[first_index:index])
                        first = True
                    index += 1
                    
    if len(items[-1]) == 0:
        items.pop(-1)           
    for pulsar in items:
        info = [param[0] for param in pulsar]
        if "ELONG" not in info and "ELAT" not in info:
            try:
                RAJ_index = info.index("RAJ")
            except:
                RAJ_index = info.index("RA")
            try:
                DECJ_index = info.index("DECJ")
            except:
                DECJ_index = info.index("DEC")
            RAJ = pulsar[RAJ_index][1]#.split(":")
            DECJ = pulsar[DECJ_index][1]#.split(":")
            coords = coordfix.SkyCoord(RAJ , DECJ, unit=(u.hourangle, u.deg), frame = "fk5")
            ELONG = coords.ra.degree
            ELAT = coords.dec.degree
            pulsar.append(['ELONG', ELONG])
            pulsar.append(['ELAT', ELAT])
    
    return items

def make_small_catalog(data, coords, catalog, **kwargs):
    """
    Makes a catalog of the obejects at the given coordinates. 
    Removes items that are not listed in the catalog
    (assuming they are not of interest)
    """
    properties = []
    new_data = []
    new_coords = []
    if kwargs["mode"] == "bands" or kwargs["mode"] == 'sliced':
        for i in range(len(data)):
            new_data.append([])
    number_not_found = 0
    number_found = 0
    for coord in coords: 
        found = False
        for item in catalog:
            ELONG_cat = float(item[[param[0] for param in item].index("ELONG")][1])
            ELAT_cat = float(item[[param[0] for param in item].index("ELAT")][1])
            ELONG_obj = coord[0]
            ELONG_err =  coord[1] + 0.075
            ELAT_obj = coord[2]
            ELAT_err =  coord[3] + 0.075
            # switch to 0, 360 degree mode (instead of -180, 180)
            if ELONG_obj < 0:
                ELONG_obj += 360

            if  ELONG_obj - ELONG_err < ELONG_cat < ELONG_obj + ELONG_err \
                and ELAT_obj - ELAT_err < ELAT_cat  < ELAT_obj + ELAT_err:
                properties.append(item)
                if kwargs["mode"] == "average":
                    new_data.append(data[coords.index(coord)])
                else:
                    for i in range(len(data)):
                        new_data[i].append(data[i][coords.index(coord)])
                new_coords.append(coord)
                found = True
                number_found += 1
                break

        if not found:
            number_not_found += 1
            if number_not_found % 1000 == 0:
                print number_not_found, " objects are not found"
            #print "No match found in catalog at coordinates: ", coord
    
    small_catalog = properties
    global small_catalog 
    
                
    print number_found, " objects matched the catalog"
    print number_not_found, " objects did not match the catalog"
    
    return new_data, new_coords, properties

def find_property(catalog, *args):
    """
    Finds properties of all objects. Poperties to be found listed in args.
    Returns a list of lists containing the properties if more than one argument 
    is provided. Otherwise a single list is returned.  
    """
    information = []
    for item in catalog:
        if len(args) > 1:
            information.append([])
            for arg in args:
                try:
                    information[-1].append(item[[i[0] for i in item].index(arg)][1])
                except:
                    information[-1].append(None)
                    #print arg, "is not a known property of this object"
        elif len(args) == 1:
            try:
                information.append(item[[i[0] for i in item].index(args[0])][1])
            except:
                information.append(None)
                print args[0], "is not a known property of this object"
                
    return information 

def save_catalog(catalog, savename):
    """
    Saves the small catalog for later usage. Can be used to speed up the process
    """
    catalog_file = open(savename, 'w')
    for pulsar in catalog:
        print "hoi"
        for item in pulsar:
            line = ""
            number = 0
            for line_item in item:
                number += 1
                line += str(line_item)
                while len(line) < 20 * number:
                    line += " "
            catalog_file.write(line + "\n")
        catalog_file.write("@" + "-" * 40 + "\n")
    catalog_file.close()
    
    print savename, "saved" 
       
def cat_hist_plot(catalog, *args):
    """
    looks up the catalog values for the flux density of the pulsars within the 
    detection range. 
    """
    
    for arg in args:
        the_list = []
        values = find_property(catalog, arg)
        for i in range(len(values)):
            if values[i]!= None:
                the_list.append(float(values[i]))
        plt.figure(arg)
        plt.title(arg)
        plt.hist(the_list, bins = 30)
        plt.show()
                        
### Determine Dates ------------------------------------------------------------

def determine_day (filename):
    """
    Determines how many days it has been since the first observation
    """
    #day = 0
    if filename[6:11] == "Feb10" or filename[4:9] == "Feb10":
        day = 0
    elif filename[6:11] == "Mar10" or filename[4:9] == "Mar10":
        day = 28
    elif filename[6:11] == "Mar24":
        day = 28 + 14
    elif filename[6:11] == "Apr21":
        day = 28 + 14 + 31 - 3
    elif filename[6:11] == "May19":
        day = 28 + 14 + 28 + 28
    elif filename[6:11] == "Jul13":
        day = 28 + 14 + 28 + 28 + 25 + 30
    elif filename[6:11] == "Jan15" or filename[4:9] == "Jan15":
        day = 365 - 31 + 5
    elif "run1" in filename:
        day = 0
    elif "run2" in filename:
        day = 28
    elif "run3" in filename:
        day = 28 + 14
    elif "run4" in filename:
        day = 28 + 14 + 28
    elif "run5" in filename:
        day = 28 + 28 + 28 + 14
    elif "run6" in filename:
        day = 28 + 28 + 28 + 30 + 25
    elif "run7" in filename:
        day = 181
    elif "run8" in filename:
        day = 181 + 79
    elif "run9" in filename: 
        day = 339
    elif "July14" in filename:
        day = 500
    elif 'PT' in filename:
        day = 0
    else:
        try: day = int(filename[32:35])
        except: 
            try: day = int(filename[32:34])
            except: day = -10
    return day

def julian_date(date):
    """
    Takes input in the form "dd-mm-yyyy"
    """
    date = date.split("-")
    day = int(date[0])
    month = int(date[1])
    year = int(date[2])
    
    a = (14 - month) / 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    
    MJD = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
    
    return MJD    

def day2sec(data, slicesize):
    """
    Changes the "day" in sliced data to the correct time in seconds
    """
    for bandno in range(len(data)):
        for pulsarno in range(len(data[bandno])):
            seconds = []
            for part in data[bandno][pulsarno][2]:
                seconds.append(part * slicesize)
            data[bandno][pulsarno][2] = seconds
    return data
### Collect Data ---------------------------------------------------------------  
         
def general(folders, name, min_flux):
    """
    Collects the fluxes from the csv files concerning the averaged images
    """  
    
    tussenstap = "D:\Project data\ "
    pulsars = []
    coords = []

    if min_flux == 0: 
        no_min_flux = True
    else:
        no_min_flux = False
    
    for folder in folders:
        print "Now working in ", folder
        os.chdir(tussenstap[:-1] + folder + tussenstap[2])
        filenames = glob.glob(os.path.expanduser(name))

        for filename in filenames:            
            data = open(filename)
            
            for line in data:
                line = line.split(", ")
                flux = line[-4]
                ferr = line[-3]
                day = determine_day(filename)
                # Skip header line
                if flux != "int_flux" and (abs(float(flux)) > min_flux or no_min_flux):
                    coord = (round(float(line[0]), 3),
                             round(float(line[1]), 3),
                             round(float(line[2]), 3),
                             round(float(line[3]), 3))
                    num = 0
                    new = True
                    for old_coord in coords:
                        # coordinates returned by pyse.py may vary
                        # 3 decimals proved to be too precise (saving one 
                        # object as multiple objects), but 1 decimal might be 
                        # not precise enough (might make one object of close 
                        # objects). 
                        if round(coord[0], 1) == round(old_coord[0], 1) and \
                           round(coord[2], 1) == round(old_coord[2], 1):
                            pulsars[num][0].append(float(flux))
                            pulsars[num][1].append(float(ferr))
                            pulsars[num][2].append(day)
                            new = False
                            break
                        num += 1
                    
                    if new: 
                        coords.append(coord)
                        pulsars.append([[float(flux)],[float(ferr)],[day]])


    print len(pulsars), "objects were found in", name
    os.chdir(tussenstap[:-1])
    return pulsars, coords

### Plotting -------------------------------------------------------------------

def band_plotting (data):
    """
    Plots all the pulsars in seperate figures and in separate bands
    """

    names = find_property(small_catalog, 'PSRJ', 'PSRB')
    bandfreq = [124, 149, 156, 185]
    index = 0
    for name in names:

        plt.figure(name[0])
        bandno = 0
        for band in data:
            plt.errorbar(band[index][2],
                     band[index][0],
                     yerr = band[index][1],
                     fmt = 'o',
                     label = "Freq: " + str(bandfreq[bandno]) + "MHz")
            plt.legend()
            bandno += 1
        
        plt.title("Flux variation for different frequency bands")
        plt.xlim(-10, 350)

        plt.show()
        index += 1
                      
def average_plot(data):
    """
    Plots the data from averaged fits files
    """

    names = find_property(small_catalog, 'PSRJ', 'PSRB')
    for name in names:
        number = names.index(name)
        try:
            plt.figure(name[0])
            #if name[0] != None:
            #    plt.title("Flux density variations of " + name[0])
            #elif name[1] != None:
            #    plt.title("Flux density variations of " + name[1])
            plt.errorbar(data[number][2],
                        data[number][0],
                        yerr = data[number][1],
                        fmt = 'o')
            plt.xlim(data[number][2][0] - 5)
            plt.xlabel("Time (Days)")
            plt.ylabel("Flux density (Jy)")
        except:
            continue
    plt.show()
        
def spectrum_plot(spec_data, *args):
    """
    plot the spectra from the modified data set
    """
    names = find_property(small_catalog, 'PSRJ', 'PSRB')
    for name in names:
        nummer = names.index(name)
        if name[0] != None:
            plt.figure(name[0])
        #    plt.title("Spectra of " + name[0])
        else:
            plt.figure(name[1])
        #    plt.title("Spectra of " + name[1])
        for day_data in spec_data[nummer]:
            if 'label' in args:
                plt.errorbar(day_data[1],
                         day_data[2],
                         yerr = day_data[3],
                         label = "Day: " + str(day_data[0]),
                         fmt = 'o')
            else: 
                plt.errorbar(day_data[1],
                         day_data[2],
                         yerr = day_data[3],
                         fmt = 'o')
            plt.legend()
            plt.xlabel("Frequency (MHz)")
            plt.ylabel("Flux density (Jy)")
            plt.show()
            
def spectrum_fit_plot(spec_data, fit_params, fit_function, *args, **kwargs):
    """
    Plots the spectra from the modified data set
    """
    names = find_property(small_catalog, 'PSRJ', 'PSRB')
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] * 3
    for name in names:
        nummer = names.index(name)
        plt.figure(name[0])
        #plt.title("Spectrum of the pulsar with powerlaw fit")
        day_no = 0
        for day_data in spec_data[nummer]:
            xdata = np.arange(120, 190, 1)
            ydata = fit_function(xdata, *fit_params[nummer][day_no])
            if 'yscale'in kwargs:
                plt.yscale(kwargs['yscale'])
            if 'xscale'in kwargs:
                plt.xscale(kwargs['xscale'])
            plt.grid(True)
            if 'label' in args:
                plt.plot(xdata, ydata, 
                         label = "Fit day: " + str(day_data[0]),
                         color = colors[day_no])
                plt.errorbar(day_data[1],
                         day_data[2],
                         yerr = day_data[3],
                         label = "Day: " + str(day_data[0]),
                         fmt = 'o',
                         color = colors[day_no])
            else: 
                plt.plot(xdata, ydata,  
                         color = colors[day_no])
                plt.errorbar(day_data[1],
                         day_data[2],
                         yerr = day_data[3],
                         fmt = 'o',
                         color = colors[day_no])
            plt.legend()
            plt.xlabel("Frequency (MHz)")
            plt.ylabel("Flux density (Jy)")
            plt.show()
            day_no += 1
            
def lots_of_plotting(data, cat, *args, **kwargs):
    """
    Plots the data in several subplots in one figure.
    The data and catalog need to be in the same order (this script does that)
    It will make multiple figures if not all plots fit on one figure. 
    Only plots 'band' data in spectrum format
    args :
        'save' Saves the figures as "Lots_of_plots_'mode'_x_to_y.png" with x the
        first plotted item (OR NUMBER OF THE PLOT/INDEX OF THE OBJECT IN DATA 
        AND THE CATALOG) and y the last. 'mode'is the mode used. 
    kwargs:
        rows: Integer, The number of rows, default = 5
        cols: Integer, The number of collumns, default = 4
        mode: String, "average" or "bands"
        figsize: tuple of floats or integers (width, height). default = (10, 14)
                 
    """
    names = find_property(cat, 'PSRJ', 'PSRB')
    try: rows = kwargs["rows"]
    except: 
        print "No number of rows defined, using 5"
        rows = 5
    try:
        cols = kwargs["cols"]
    except:
        print "No number of columns defined, using 4"
        cols = 4
    ## For consistent coloring of different epcohs.     
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] * 5
    
    try: 
        mode = kwargs["mode"]
        if mode != 'bands' and mode != 'average':
            print "'", mode, "'", "is an invalid mode, using 'average'"
            mode = 'average'
    except: 
        print "No mode defined, using 'average'"
        mode = "average"
    try: size = kwargs["figsize"]
    except: size = (10, 14)
    
    if cols * rows >= len(data):
        f, axarr = plt.subplots(rows, 
                                cols, 
                                figsize = size)
        for i in range(len(data)):
            plot_row = i / cols
            plot_col = i % cols
            name = names[i]
            
            if name[1] != None:
                axarr[plot_row, plot_col].set_title(name[1],fontsize = 11)
            elif name[0] != None:
                axarr[plot_row, plot_col].set_title(name[0],fontsize = 11)
                
            if mode == 'average':
                x_min = min(data[i][2]) - 10
                x_max = max(data[i][2]) + 10
                axarr[plot_row, plot_col].locator_params(axis='x',nbins=4)
                axarr[plot_row, plot_col].set_xlim([x_min, x_max])
                axarr[plot_row, plot_col].errorbar(data[i][2],
                                                   data[i][0],
                                                   yerr=data[i][1],
                                                   fmt = 'o')
            elif mode == 'bands':
                axarr[plot_row, plot_col].set_xlim([100, 200])
                plt.sca(axarr[plot_row, plot_col])
                plt.xticks([100, 150, 200])
                color_num = 0
                for day in data[i]:
                    axarr[plot_row, plot_col].errorbar(day[1],
                                                       day[2],
                                                       yerr = day[3],
                                                       fmt = "o",
                                                       color = colors[color_num])
                    color_num += 1
                    
        f.tight_layout()
        if 'save' in args:
            f.savefig("Lots_of_plots_" + mode + "_1_" + str(len(data)) + ".png") 
                                                                          
    elif cols * rows < len(data):
        for fnum in range(len(data) / (rows * cols) + 1):
            pnum_0 = fnum * cols * rows
            pnum_1 = (fnum + 1) * cols * rows
            f, axarr = plt.subplots(rows, 
                                    cols, 
                                    figsize = size)
            
            if pnum_1 > len(data):
                pnum_1 = len(data)
            for i in range(pnum_0, pnum_1):
                plot_row = (i - fnum * cols * rows) / cols
                plot_col = i % cols
                name = names[i]
                if name[1] != None:
                    axarr[plot_row, plot_col].set_title(name[1],fontsize = 11)
                elif name[0] != None:
                    axarr[plot_row, plot_col].set_title(name[0],fontsize = 11)
                
                if mode == 'average':
                    x_min = min(data[i][2]) - 10
                    x_max = max(data[i][2]) + 10
                    axarr[plot_row, plot_col].locator_params(axis='x',nbins=4)
                    axarr[plot_row, plot_col].set_xlim([x_min, x_max])
                    axarr[plot_row, plot_col].errorbar(data[i][2],
                                                            data[i][0],
                                                            yerr=data[i][1],
                                                            fmt = 'o')
                elif mode == 'bands':
                    axarr[plot_row, plot_col].set_xlim([100, 200])
                    plt.sca(axarr[plot_row, plot_col])
                    plt.xticks([100, 150, 200])
                    color_num = 0
                    for day in data[i]:
                        axarr[plot_row, plot_col].errorbar(day[1],
                                                                day[2],
                                                                yerr = day[3],
                                                                fmt = "o",
                                                                color = colors[color_num])
                        color_num += 1
                        
            f.tight_layout()
            if 'save' in args:
                f.savefig("Lots_of_plots_" + mode + "_" + str(pnum_0 + 1) + "_" + str(pnum_1) + ".png")
    
    plt.show()
                 

def sliced_graph(data):
    """
    Plots the sliced datapoints in one figure
    """
    names = find_property(small_catalog, 'PSRJ', 'PSRB')
    for name in names:
        number = names.index(name)
        for band in data:
            plt.figure(name[0])
            plt.errorbar(band[number][2],
                        band[number][0],
                        yerr = band[number][1],
                        fmt = 'o')
            plt.xlabel("Time (Seconds)")
            plt.ylabel("Flux density (Jy)")

    plt.show()

def sliced_plot(data, freqs, *args):
    """
    not in use at the moment, replaced by dynamic_spectrum_plot
    """
    
    plt.figure("Colormap2")
       
    times = data[0][0][2]
    
    
    matrix = [x[0][0] for x in data] 
    #matrix[2][2:-1] = [np.nan] * len(matrix[2][2:-1])
    matrix = rotate_matrix(matrix)
    #Inverting x axis.
    matrix = [x[::-1] for x in matrix]
    try:
        matrix, times = add_empty_space(matrix, times)
    except:
        pass
    masked_array = np.ma.array (matrix, mask=np.isnan(matrix))
    if "imshow" in args:
        cmap = plt.cm.cool
        cmap.set_bad('w',1.)
        #plt.imshow(masked_array, interpolation='none', cmap=cmap)
    if "pcolor" in args:
        times.append(2 * times[-1] - times[-2])
        freqs.append(2 * freqs[-1] - freqs[-2])
        plt.pcolor(freqs, times, masked_array, cmap='cool')
    if "contour" in args:
        plt.contour(freqs, times, masked_array)
    plt.xlabel("(Frequency subbands for now)")
    plt.ylabel("Time (s)")
    plt.colorbar()
    plt.show()

def dynamic_spectrum(data, **kwargs):
    """
    Plots the dynamic spectrum
    """
    ## Check frequency length and data length
    if 'scale' in kwargs:
        scale = kwargs['scale']
    else: scale = [0, 6.5]
    freqs = kwargs['freqs']
    freq_len = 0
    for freq_list in freqs:
          freq_len += len(freq_list)
    if len(data) != freq_len:
        print "The frequencies and data do not agree (lengths are not the same)"
    
    titles = ['124 MHz', '149 MHz', '156 MHz', '185 MHz']
    big_matrix = [x[0][0] for x in data]
    times = [data[0][0][2]] * len(freqs)
    matrices = []
    
    count = 0
    for freq_list in freqs:
        matrices.append(big_matrix[count:count + len(freq_list)])
        count += len(freq_list)
        
    for i in range(len(matrices)):
        matrices[i] = rotate_matrix(matrices[i])
        matrices[i] = [x[::-1] for x in matrices[i]]
        try:
            matrices[i], times[i] = add_empty_space(matrices[i], times[i])
        except:
            pass
        matrices[i] = np.ma.array (matrices[i], mask=np.isnan(matrices[i]))
    
    for i in range(len(freqs)):
        times[i] += [2 * times[i][-1] - times[i][-2]]
        freqs[i] += [2 * freqs[i][-1] - freqs[i][-2]]

    f, axarr = plt.subplots(1, len(freqs), sharey = True)
    for i in range(len(matrices)):
        axarr[i].set_title(titles[i],fontsize = 11)
        axarr[i].set_ylim([0, 1540])
        axarr[i].set_xlim([freqs[i][0], freqs[i][-1]])
        if scale != None:
            im = axarr[i].pcolor(freqs[i], times[i], matrices[i], cmap = 'cool', vmin = scale[0], vmax = scale[1])
        else:
            im = axarr[i].pcolor(freqs[i], times[i], matrices[i], cmap = 'cool')
    axarr[0].set_ylabel("Time (seconds)")
    f.text(0.47, 0.04, 'Frequency (MHz)', ha='center')
    #f.colorbar(im, ax=axarr.ravel().tolist())        
    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.82, 0.1, 0.03, 0.8])
    f.colorbar(im, cax=cbar_ax)
    f.subplots_adjust(wspace=0)
    plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
    
    plt.show()
     
    
    
    
    
### Data cleaning and manipulation ---------------------------------------------

def rotate_matrix(a):
    """
    rotates a given matrix, *a*,  90 degrees clockwise
    """
    c = [[a[y][x] for y in range(len(a))[::-1]] for x in range(len(a[0]))]
    return c

def add_empty_space(matrix, times):
    """
    Adds the moment of calibrations
    Adds the absent data to the matrix as NaN values
    Assumes time starts at 0 and all steps are the same size. 
    Assumes only one empty area
    """
    new_times = range(0, times[-1], times[1])
    to_add = (len(new_times) - len(times)) * [[np.nan] * len(matrix[0])]
    
    for i in range(len(times) - 1):
        if times[i + 1] - times[i] != times[1]:
            sep_index = i + 1
            break
    new_matrix = matrix[:sep_index] + to_add + matrix[sep_index:]
    
    return new_matrix, new_times

def add_missing_data(data):
    """
    Makes sure all data is the same lenght, adds NaN as data for "missing" data
    Also assumes equal sized steps and first time step at 0. 
    """
    length = max([len(band[0][0]) for band in data])
    for band in range(len(data)):
        if len(data[band][0][0]) == length:
            correct_times = data[band][0][2]
            break
    for band in range(len(data)):
        if len(data[band][0][0]) != length:
            data[band][0][0] += [np.nan] * (length - len(data[band][0][0]))
            data[band][0][1] += [np.nan] * (length - len(data[band][0][1]))
            data[band][0][2] = correct_times
        
    return data


def spectrum_data(data):
    """
    changes the data format to frequency dependent fluxes, for separate days
    """

    names = find_property(small_catalog, 'PSRJ', 'PSRB')
    spec_data = []
    for name in names:
        nummer = names.index(name)
        spec_data.append([])
        count = 0
        for day in data[0][nummer][2]:
            spec_data[nummer].append(
                     [day, [124, 149, 156, 185],
                     [data[0][nummer][0][count],
                      data[1][nummer][0][count],
                      data[2][nummer][0][count],
                      data[3][nummer][0][count]],
                      
                     [data[0][nummer][1][count],
                      data[1][nummer][1][count],
                      data[2][nummer][1][count],
                      data[3][nummer][1][count]]])
            count += 1

    return spec_data
    
def clean_avg_data (data):
    """
    Removes the data points very close to 0 (-0.001 < data point < 0.001), because
    the images contain empty space outside the 
    """
    for pulsar in data:
        nummer = data.index(pulsar)
        index_to_remove = []
        count = 0
        for flux in pulsar[0]:
            
            if abs(float(flux)) < 0.001:
                index_to_remove.append(count)
            count += 1
        removed = 0
        for i in index_to_remove:
            data[nummer][0].pop(i - removed)
            data[nummer][1].pop(i - removed)
            data[nummer][2].pop(i - removed) 
            removed += 1   
    return data

def clean_double_data (data):
    """
    Looks at the data and removed the data that is a result of overlapping
    images. Only saves the data of the highest quality. No functional data is
    lost. Cleans data for average mode.
    """
    ## note: The data that is saved will be from one image. 
    ## The data is first sorted to make sure of that. 
    
    for pulsar in data:
        # list the epochs and their index in the list
        days = [[]]
        indexes = [[]]
        for i in range(len(pulsar[0])):
            found = True
            j = 0
            while found:
                try:
                    if pulsar[2][i] not in days[j]:
                        days[j].append(pulsar[2][i])
                        indexes[j].append(i)
                        found = False #It was not in  this list
                    else: j += 1
                except IndexError: #create a new list when the last list is passed
                    days.append([pulsar[2][i]])
                    indexes.append([i])
                    found = False #It is not in this new list
        if len(days) > 1:
            mean_sig = []
            for index_list in indexes:
                flux_mean = np.mean([pulsar[0][n] for n in index_list])
                error_mean = np.mean([pulsar[1][n] for n in index_list])
                mean_sig.append(flux_mean / error_mean)
            # The information to save
            save_list = mean_sig.index(max(mean_sig)) #index in the list of indexes
            index_list = indexes[save_list] # list of indexes
            #change the data
            data[data.index(pulsar)] = [[pulsar[0][n] for n in index_list],
                                        [pulsar[1][n] for n in index_list],
                                        days[save_list]]
    
    return data

def remove_empty_data(data, catalog):
    """
    Removes the objects that have lost all data in the cleaning process. 
    These are pulsars that lie within an image, but outsude the data containing
    part of the image.
    """
    remove = []
    for i in range(len(data)):
        if len(data[i][0]) < 1:
            remove.append(i)
    removed = 0
    for i in remove:
        data.pop(i - removed)
        print catalog[i][0][0], ": ", catalog[i][0][1], "is being removed due to empty data lists."
        catalog.pop(i - removed)
        
        removed += 1
    
    return data, catalog                

def making_arrays(data):
    """
    changes all lists to numpy arrays, works for both the bands and the spectrum
    data formats. Also works for the average value data.
    """
    num = 0
    #Spectrum
    try:
        if type(data[0][0][0]) == int:
            for pulsar in data:
                data[num][0][1] = np.array(data[num][0][1])
                data[num][0][2] = np.array(data[num][0][2])
                data[num][0][3] = np.array(data[num][0][3])
                num += 1
    except: 
        pass
    #Band 
    try:      
        if type(data[0][0][0]) == list: 
            for pulsar in data:
                data[num][0] = np.array(data[num][0])
                data[num][1] = np.array(data[num][1])
                data[num][2] = np.array(data[num][2])
                num += 1
    except: 
        pass
    #Average
    try:        
        if type(data[0][0][0]) == float:
            for pulsar in data:
                data[num] = np.array(data[num])
                num += 1
    except: 
        pass
        
    return data
                            
def spec_clean(data):
    """
    Selects the best data from the different images, discards the other data
    """
    clean_data = []
    
    for pulsar in data:
        clean_data.append([])
        count1 = 0 # keeps track of the index of day1
        double_days = []
        for day1 in pulsar:
            double_days.append([])
            count2 = 0 # keeps track of the index of day2
            for day2 in pulsar:
                # if the days(dates) are equal, but the data is not
                if day1[0] == day2[0] and not day1 == day2:
                    if count1 not in double_days[-1]:
                        double_days[-1].append(count1)
                    if count2 not in double_days:
                        double_days[-1].append(count2)
                count2 += 1
            if len(double_days[-1]) == 0:
                double_days[-1].append(count1)
            count1 += 1    
        for i in range(len(double_days)):
            double_days[i] = list(np.sort(double_days[i]))
        for i in double_days:
            if double_days.count(i) > 1:
                double_days.pop(double_days.index(i))
        for indices in double_days:
            values = []
            for i in indices:
                values.append(np.mean(pulsar[i][2]) / np.mean(pulsar[i][3]))

            
            clean_data[-1].append(pulsar[indices[values.index(np.max(values))]])

    return clean_data
    
### Dispersion measure and variation -------------------------------------------

def find_DM():
    
    DM_values = find_property(small_catalog, 'DM')
    DM_list = []
    
    for DM in DM_values:
        if DM != None:
            DM_list.append(float(DM))
        else: DM_list.append(None)

    return DM_list
                
def find_var(data):
    """
    calculates the variance of the flux of the pulsars. Returns a list
    """
    var_list = []
    for pulsar in data:
        var_list.append(np.var(pulsar[0]))  
    return var_list

def find_var2(data):
    """
    calculates the relative variance of the fluxes, takes errors into account
    (kind of). Note: Data has to be in np.array format (not list)
    """
    var_list = []
    for pulsar in data:
        fmean = np.mean(pulsar[0])
        var = np.mean(abs((pulsar[0] - fmean) / pulsar[1]))
        var_list.append(var)
        
    return var_list

def find_var3(data):
    """
    calculates the relative variance of the fluxes, takes errors into account
    (kind of). Note: Data has to be in np.array format (not list)
    """
    var_list = []
    for pulsar in data:
        fmean = np.mean(pulsar[0])
        error = pulsar[1]
        flux = pulsar[0]
        var = np.mean(abs((flux - fmean)**2 / (error * fmean)))
        var_list.append(var)
        
    return var_list

def find_var4(data):
    """
    Calculates variability as the flux coeffecient variation, first option from
    http://tkp.readthedocs.io/en/latest/userref/structure/stages/transient.html
    """
    var_list = []
    for pulsar in data:
        fmean = np.mean(pulsar[0])
        flux = pulsar[0]
        N = len(flux)
        var =(1/fmean) * np.sqrt((N * (np.mean(pulsar[0]**2) - fmean**2))/(N-1))
        var_list.append(var)
    return var_list                                  

def find_var5(data):
    """
    Calculates the variability of the flux difined based on reduced chi squared
    statistics. From: 
    http://tkp.readthedocs.io/en/latest/userref/structure/stages/transient.html
    Note: This method has a flawed significance calculation 
    """
    var_list = []
    for pulsar in data:
        error = pulsar[1]
        flux  = pulsar[0]
        N     = len(flux)
        weighted_mean_flux = np.sum(flux / error) / np.sum(1. / error)
        var = 1. / (N - 1) * np.sum((flux - weighted_mean_flux)**2 / error**2) 
        var_list.append(var)
    return var_list                  
   

def print_variation (data, small_catalog, function):
    """
    Prints the variation of the pulsars in the data calculated with the 
    specified function. The specified function has to take only data as input 
    and return a list of variabilities in the same order.
    """
    var_list = function(data)
    for i in range(len(var_list)):
        print small_catalog[i][0][1], "\t", small_catalog[i][1][1]
        print "Variation: ", var_list[i]

def DM_var_plot(data):
    DM_list = find_DM()
    var_list = find_var(data)
    var_list2 = find_var2(data)
    var_list3 = find_var3(data)
    var_list4 = find_var4(data)
    var_list5 = find_var5(data)
    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex = True)
    ax1.plot(DM_list, var_list,'o', label = "Numpy variance")
    ax1.legend()
    ax2.plot(DM_list, var_list2,'o', label = "Custom variance")
    ax2.legend()
    ax3.plot(DM_list, var_list3,'o', label = "Custom variance 2")
    ax3.legend()
    ax4.plot(DM_list, var_list4,'o', label = "Coefficient of variation")
    ax4.legend()
    ax5.plot(DM_list, var_list5,'o', label = "Reduced chi^2")
    ax5.legend()
    plt.xlabel("DM (pc/cm$^{3}$)")
    plt.ylabel("Variance")
    plt.show()
    
def DM_var_plot2(data):
    DM_list = find_DM()
    var_list5 = find_var5(data)
    plt.plot(DM_list, var_list5, 'o')
    plt.xlabel("DM (pc/cm$^{3}$)")
    plt.ylabel("$\eta$")
    plt.yscale('log')
    plt.show()
### Object Selecion ------------------------------------------------------------

def detected_only(data, coords, filename, mode):
    """
    Checks the given file for objects of interst, all other objects are removed
    from the collected data and temporary (small_)catalog 
    """
    global small_catalog
    names = find_property(small_catalog, 'PSRJ', 'PSRB')

    detected_pulsars = read(filename)
    print len(detected_pulsars), "interesting objects are listed"
    new_coords = []
    new_data = []
    new_small_catalog = []
    print "Names: ", len(names)
    print "Data: ", len(data)
    if mode == 'bands' or mode == 'sliced':
        for band in range(len(data)):
            new_data.append([])
            
            for name in names:
                if name[0] in detected_pulsars:
                    new_data[band].append(data[band][names.index(name)])
                    if band == 0: #add the coordinates only once
                        new_coords.append(coords[names.index(name)])
                        new_small_catalog.append(small_catalog[names.index(name)])
        print len(new_data[0]), "interesting objects are found" 
    else:
        for name in names:
            if name[0] in detected_pulsars:
                new_data.append(data[names.index(name)])
                new_coords.append(coords[names.index(name)])
                new_small_catalog.append(small_catalog[names.index(name)])
        print len(new_data), "interesting objects are found"
    small_catalog = new_small_catalog    
    return new_data, new_coords

### Fitting --------------------------------------------------------------------

def fitting (data, fit_function, *args):
    """
    Fits the spectra to the fit_function, using curve_fit. Returns a list of 
    parameters for every pulsar on every date
    """
    params = []
    covars = []
    for pulsar in data:
        params.append([])
        covars.append([])
        for day in pulsar:
            param, pcov = curve_fit(fit_function,
                                    day[1], 
                                    day[2], 
                                    sigma=day[3], 
                                    maxfev = 5000000, 
                                    absolute_sigma = True)
                                    
            perr = np.sqrt(np.diag(pcov))
            params[-1].append(param)
            covars[-1].append(perr)
    if 'covar' in args:
        return params, covars
    else:
        return params    
            
def fit_func(x, a, b):
    return a * x**b

def power(x, a, b):
    return a * x**b

def power_test(x, a, b, c):
    return a + b * x**c

def print_fit_params(data, catalog, fit_function):
    """
    Prints the fit parameters of the fit_function on the spectrum data.
    """
    params, covars = fitting(data, fit_function, 'covar')
    
    for i in range(len(data)):
        print catalog[i][0][1], "\t", catalog[i][1][1]
        for j in range(len(data[i])):
            print "Day: %s" % data[i][j][0]
            print "a: ", params[i][j][0], "\t Standard deviation: ", covars[i][j][0]
            print "b: ", params[i][j][1], "\t Standard deviation: ", covars[i][j][1]
            if len(params[i][j]) == 3:
                print "c: ", params[i][j][2], "\t Standard deviation: ", covars[i][j][2], "\n"
        print "\n"
            
        

### Master Functions -----------------------------------------------------------
    
def run_bands(directories, filenames, big_catalog, *args, **kwargs):
    """
    Runs all seperate parts of the code if the arguments are present. Plots
    the result. Returns found values
    Kwargs are:
        min_flux = float (minimum value of the flux (detection threshold))
        
        detected = file with names of specific items of interest. Only these 
                   will remain in the data.
        
        fit = function (name of a function that can be used to fit the data)   
    """
    data = []
    if "min_flux" in kwargs:
        min_flux = kwargs["min_flux"]
    else: min_flux = 0
    
    #collecting all data from csv files
    for filename in filenames:
        band_data, coords = general(directories, filename, min_flux)
        data.append(band_data)
  
    data, coords, small_catalog = make_small_catalog(data, coords, big_catalog, mode = 'bands')
    global small_catalog
    
    if "detected" in kwargs:
        data, coords = detected_only(data, coords, kwargs['detected'], kwargs['mode'])
    # making it easier to plot a spectrum

    data = spectrum_data(data)
    data = spec_clean(data)
    data = making_arrays(data)
    #print_fit_params(data, small_catalog, power)
    if 'plotlots' in args:
        lots_of_plotting(data, small_catalog, *args, **kwargs)
           
    if 'fit' in kwargs:
        fit_function = kwargs['fit']
        params = fitting(data, fit_function)
        
        if 'SPINDX' in args:
            compare_SPINDX(params, small_catalog)
            
        if 'plot' in args:
            spectrum_fit_plot(data, params, fit_function, *args, mode = 'average')
        
        return data, coords, params
        
    else: 
        if 'plot' in args:
            spectrum_plot(data, *args)
            
        return data, coords

        
def run_avg(directories, filenames, big_catalog, *args, **kwargs):
    """
    Collects data from the averaged images
    """
    data = []
    if "min_flux" in kwargs:
        min_flux = kwargs["min_flux"]
    else: min_flux = 0
    #collecting all data from csv files
    for filename in filenames:
        d, coords = general(directories, filename, min_flux)
        data += d
        
    # Creating a smaller catalog to be used in functions    
    data, coords, small_catalog = make_small_catalog(data, coords, big_catalog, mode = 'average')
    global small_catalog
    
    if 'SPINDX' in args:
        fluxes = calc_flux(small_catalog, big_catalog)
        compare_fluxes(fluxes, data, small_catalog)
    
    
    if "detected" in kwargs:
        data, coords = detected_only(data, coords, kwargs['detected'], 'average')
    
    if 'clean' in args:
        data = clean_avg_data(data)
        data = clean_double_data(data)
        # Dubious function
        data, small_catalog = remove_empty_data(data, small_catalog)
    
    if 'DM' in args:
        data = making_arrays(data)
        DM_var_plot(data)
    if "DM2" in args:
        data = making_arrays(data)
        DM_var_plot2(data)
        
    if 'plotlots' in args:
        lots_of_plotting(data, small_catalog, *args, **kwargs)
        
    if 'plot' in args:
        average_plot(data) 
      
    return data, coords

def run_sliced(directories, filenames, big_catalog, *args, **kwargs):
    """
    Runs the sliced mode
    """
    data = []
    if "min_flux" in kwargs:
        min_flux = kwargs["min_flux"]
    else: min_flux = 0
    
    if "slicesize" in kwargs:
        slicesize = kwargs["slicesize"]
    else: slicesize = 30
    
    #collecting all data from csv files
    for filename in filenames:
        sliced_data, coords = general(directories, filename, min_flux)
        data.append(sliced_data)
  
    data, coords, small_catalog = make_small_catalog(data, coords, big_catalog, mode = 'sliced')
    if "detected" in kwargs:
        data, coords = detected_only(data, coords, kwargs['detected'], kwargs['mode'])
    
    #Below here things function properly   
    data = day2sec(data, slicesize)
    data = add_missing_data(data)
    freqs = kwargs["freqs"]
    sliced_graph(data)
    #sliced_plot(data, freqs, *args)
    dynamic_spectrum(data, **kwargs)
    
    global small_catalog
    return data, coords
        
def collect_data(directories, filenames, *args, **kwargs):
    """
    The master function. (should be all you need)
     
    Collects the data from the given csv files. The csv files have to be in the 
    format pyse.py returns 
    
    
    Returns the data in the form:

        for averaged images:
            list of pulsars[(list of fluxes per day[], list of errors per day[],
            list of days[] , pulsar name) ]
        for bands:
            list of pulsars[ (day, list of frequencies[], list of fluxes[], list
            of errors[], pulsar name) ]
        pulsar name as tuple (PRSJ, PSRB)
        
    input: 
        directories: list of strings (the directories of the csv files)
        filenames: list of strings (names of the csv files, "*partial name*"
                   can be used.)
                   If mode == bands every name will create a new band, so use 
                   "*partialname*" notation for multiple files in one band.
       args:
           plot: Plots the results
           
           clean: Cleans data to prevent overlap, saves the best data, discards
                  the other data.
           DM: Plots the variability of the flux as a function of the DM 
               (Dispersion measure)
               
           label: Puts labels in the plots (can be rather large)
           
           SPINDX: Prints calculated(expected) flux densities as a result of the
                   given spectral indices and flux densities (at a certain
                   frequency) found in the catalog. If a pulsar does not have a
                   spectral index listed the average value of all pulsars is 
                   used. Works only for "average" mode.
           plotlots: Plots the data in giant figures. The size of the figure,
                   amount of rows and columns can be defined with kwargs. 
                   See help(lots_of_plotting) for more help. 
           save:   saves the giant figures of 'plotlots' as 
                   "Lots_of_plots_x_y.png" with x the index of the first item
                   in the figure and y the last. 
                 
       kwargs: 
           catalog: string, The filename of the catalog that will be used for 
                the determination of the properties of the sources. 
                (the file has to be in the working directory or the entire path
                has to be specified)
                
           detected: string (filename of txt file with the name of pulsars of 
                interest. Only works with selected pulsars)
                
           min_flux: float, requires the data to be about a minimum flux value
           
           mode: string, options are:
                bands: saves the data in a band format to plot spectra
                average: saves data in a single band, no different frequencies
                Default is average.
           fit: function (fits the spectra with the given function, only works 
                with bands)
           xscale and yscale: string like the input for matplotlib figures.
                If not specified matplotlib default 'linear' is used. 
                Only affects the fitted figures. 
    """
    print "Reading catalog"
    big_catalog = read_database(kwargs['catalog'])
    print "Found ", len(big_catalog), " objects in catalog"
    if 'mode' in kwargs:
        if kwargs['mode'] == 'bands':
            if 'fit' in kwargs:
                data, coords, params = run_bands(directories, 
                                                 filenames, 
                                                 big_catalog,
                                                 *args, 
                                                 **kwargs)
                return data, params, small_catalog
                                                 
            else: data, coords = run_bands(directories,
                                           filenames, 
                                           big_catalog, 
                                           *args, 
                                           **kwargs)
            
        elif kwargs['mode'] == 'average':
            data, coords = run_avg(directories,
                                   filenames,
                                   big_catalog,
                                   *args,
                                   **kwargs)
        elif kwargs['mode'] == 'sliced':
            data, coords = run_sliced(directories, 
                                      filenames,
                                      big_catalog,
                                      *args,
                                      **kwargs)
        else: 
            print kwargs['mode'], "is not a possible mode, using 'average'."
    else:
        data, coords = run_avg(directories,
                               filenames,
                               big_catalog, 
                               *args,
                               **kwargs)
                
    return data, small_catalog

### Spectral index -------------------------------------------------------------

def S150 (spindx, freq, S_freq):
    """
    Calculates the intensity (flux density) at 150 MHz using S ~ F^spindx
    """
    return (150. / freq)**spindx * S_freq

def average_SPINDX (catalog):
    """
    Calculates the average Spectral index of the pulsars in the catalog
    """
    SPINDX_list = []
    
    for pulsar in catalog:
        for entry in pulsar:
            if entry[0] == "SPINDX":
                SPINDX_list.append(float(entry[1]))
    
    print len(SPINDX_list), "Spectral indices found"
    print np.mean(SPINDX_list), "is the mean value"
    print np.std(SPINDX_list), "is the standard deviation"
    
    return np.mean(SPINDX_list), np.std(SPINDX_list)
    
def calc_flux (small_catalog, big_catalog):
    """
    Calculates the expected flux at 150 MHz with flux densities and spectral 
    index from the catalog. For pulsars without an spectral index listed the 
    average spectral index of all pulsars is used. 
    """
    backup_SPINDX = average_SPINDX(big_catalog)
    
    calced_fluxes = []
    for pulsar in small_catalog:
        entries = [i[0] for i in pulsar]
        calced_fluxes.append([])
        if "SPINDX" in entries:
            SPINDX = float(pulsar[entries.index("SPINDX")][1])
        else:
            SPINDX = backup_SPINDX[0]
        calced_fluxes[-1].append("SPINDX: " + str(SPINDX))
        temp = [] # calculated fluxes
        for i in entries:
            try:
                if i[0] == "S":
                    freq = float(i[1:])
                    flux = float(pulsar[entries.index(i)][1]) / 1000
                    temp.append(S150(SPINDX, freq, flux))
            except: continue
        calced_fluxes[-1].append(temp)
                
    return calced_fluxes
    
def compare_fluxes (calced_fluxes, data, small_catalog):
    """
    Prints the measured flux(es) and the calculated flux(es). 
    """
    count = 0
    for i in range(len(calced_fluxes)):
        flux = max(data[i][0])
        erro = data[i][1][data[i][0].index(max(data[i][0]))]
        try:
            sigma = flux / erro
        except: sigma = 0
        #if sign > 5.5:
        try:
            max_calced = max(calced_fluxes[i][1])
        except ValueError: 
            max_calced = 0
            print small_catalog[i][0][0], "\t", small_catalog[i][0][1], \
                   "Had no flux densities listed in the catalog \n"
        if max_calced > 0.15:           
            print small_catalog[i][0][0], "\t", small_catalog[i][0][1]
            print small_catalog[i][1][0], "\t", small_catalog[i][1][1]
            print "Found:\t", "Flux: ", flux, "\t", "Error: ", erro, "\t", 'sigma:', sigma
            print "Calced:\t", "Mean Flux: ", np.mean(calced_fluxes[i][1]), "\t", calced_fluxes[i][0]
            print "\t", "Individual values: ", calced_fluxes[i][1], "\n"
            count += 1
    print count, "Pulsars should have (should) been detected"

def compare_SPINDX (fit_params, small_catalog):
    """
    Compares the calculated values of the spectral index with the ones in the 
    catalog
    """
    calced_indexes = []
    for i in range(len(fit_params)):
        calced_indexes.append([index[1] for index in fit_params[i]])
    cat_indexes = find_property(small_catalog, 'SPINDX')
    print len(fit_params), len(calced_indexes), len(cat_indexes)
    for i in range(len(fit_params)):
        print small_catalog[i][0][0], "\t", small_catalog[i][0][1]
        print small_catalog[i][1][0], "\t", small_catalog[i][1][1]
        print "Calculated: ", calced_indexes[i]
        print "Mean Calculated: ", np.mean(calced_indexes[i])
        print "Catalog: ", cat_indexes[i], "\n"
        
    
        
        
### Misc -----------------------------------------------------------------------
def extra_clean():
    """
    temporary function to compare results
    
    collect_data(['PT19_SAP005'],
                 ['*BAND00*', '*BAND01*', '*BAND02*', '*BAND03*'],
                 'plot',
                 'clean',
                 'label', 
                 mode = 'bands',
                 fit = power,
                 catalog = 'psrcat.db')
    """
    os.chdir("D:\Project data")
    collect_data(['PT19_SAP005_v2'],
                 ['*BAND00*', '*BAND01*', '*BAND02*', '*BAND03*'],
                 'plot',
                 'clean', 
                 mode = 'bands',
                 fit = power,
                 catalog = 'psrcat.db')
                 
def sig_test (data):
    """
    Tests the significance for every data point
    """ 
    for band in data:
        for pulsar in band:
            for i in range(len(pulsar[0])):
                print "Sig: %f" % (pulsar[0][i] / pulsar[1][i])
### End ------------------------------------------------------------------------

os.chdir("D:\Project data")
mosaics = glob.glob(os.path.expanduser("csv*"))
sig3 = glob.glob("*sig3")
sig6 = glob.glob("*6sig*")     
beams = ['Pointings_CSV_v2']
extra = ['extra_data']
bands = ['*BAND00*', '*BAND01*', '*BAND02*', '*BAND03*']
averages = ['*avg.csv']
subbands = ['*BAND00*subband0*', 
            '*BAND00*subband1*', 
            '*BAND00*subband2*', 
            '*BAND00*subband3*',
            '*BAND01*subband0*', 
            '*BAND01*subband1*', 
            '*BAND01*subband2*', 
            '*BAND01*subband3*',
            '*BAND02*subband0*', 
            '*BAND02*subband1*', 
            '*BAND02*subband2*', 
            '*BAND02*subband3*',
            '*BAND03*subband0*', 
            '*BAND03*subband1*', 
            '*BAND03*subband2*', 
            '*BAND03*subband3*']
            
subbands2 = ['*BAND00*subband0*', 
            '*BAND00*subband1*', 
            '*BAND01*subband0*', 
            '*BAND01*subband1*', 
            '*BAND02*subband0*', 
            '*BAND02*subband1*', 
            '*BAND03*subband0*', 
            '*BAND03*subband1*']
            
subbands10 = ["*BAND00*subband" + str(i) + "*" for i in range(10)]
subbands10 += ["*BAND01*subband" + str(i) + "*" for i in range(10)]
subbands10 += ["*BAND02*subband" + str(i) + "*" for i in range(10)]
subbands10 += ["*BAND03*subband" + str(i) + "*" for i in range(10)]

frqs = [[123045349.121, 123338317.871, 123631286.621, 123924255.371, 124217224.121,
          124412536.621, 124510192.871, 124607849.121, 124705505.371, 124803161.621],
         [148045349.121, 148338317.871, 148631286.621, 148924255.371, 149217224.121,
          149412536.621, 149510192.871, 149607849.121, 149705505.371, 149803161.621],
         [155076599.121, 155369567.871, 155662536.621, 155955505.371, 156248474.121,
          156443786.621, 156541442.871, 156639099.121, 156736755.371, 156834411.621],
         [183982849.121, 184275817.871, 184568786.621, 184861755.371, 185154724.121, 
          185350036.621, 185447692.871, 185545349.121, 185643005.371, 185740661.621]]
          
#frqs = np.array(frqs) / 1000000.
frqs = [[x / 1000000 for x in frq] for frq in frqs]        

### Collecting data
collect_data(mosaics, averages, 'plot', 'label', 'clean', detected = "Detected v2.txt", mode = 'average', catalog = 'psrcat.db')
#collect_data(beams, bands, 'plot', 'clean', 'label', 'DM', detected = 'Detected v2.txt', mode = 'bands', fit = power, catalog = 'psrcat.db')
#collect_data(mosaics, bands, 'plot', 'clean', 'label', detected = 'Detected v2.txt', mode = 'bands', fit = power, catalog = 'psrcat.db')

### Sliced data
#collect_data(['B0329+54_only'], subbands,'pcolor', 'plot', 'label', detected = "B0329+54.txt", mode = 'sliced', catalog = 'Small catalog averages.txt', freqs = range(16))
#collect_data(['15sec'], subbands2, 'plot','pcolor', 'label', detected = "B0329+54.txt", mode = 'sliced', catalog = 'Small catalog averages.txt', freqs = range(8))
#collect_data(['20sec'], subbands2, 'plot', 'label','contour', detected = "B0329+54.txt", mode = 'sliced', catalog = 'Small catalog averages.txt', freqs = range(8), slicesize = 20)
#collect_data(['20sec_4subbands'], subbands,'pcolor', 'plot', 'label', detected = "B0329+54.txt", mode = 'sliced', catalog = 'Small catalog averages.txt', freqs = range(16), slicesize = 20)
#collect_data(['BAND03_20sec_10subbands'], subbands10,'pcolor', 'plot', 'label', detected = "B0329+54.txt", mode = 'sliced', catalog = 'Small catalog averages.txt', freqs = range(10), slicesize = 20)
#Not in use anymore use VV 

#collect_data(['Almost_all_10subbands'], subbands10,'pcolor', 'plot', 'label', detected = "B0329+54.txt", mode = 'sliced', catalog = 'Small catalog averages.txt', freqs = frqs, slicesize = 20)
#collect_data(['Stable_source'], subbands10, mode = 'sliced', catalog = 'Fake_cat.txt', freqs = frqs, slicesize = 20)
collect_data(['Noise'], subbands10, mode = 'sliced', catalog = 'Fake_cat.txt', freqs = frqs, slicesize = 20, scale = None)

### Check detected and should be detected
# data, cat = collect_data(["Stacked"], ["*.csv"], 'SPINDX', mode = 'average', catalog = 'psrcat.db')

### Variation commands
#print_variation(data, small_catalog, find_var5)

### Spectral index commands
#run_avg(mosaics, averages, 'plot', 'clean', 'DM', detected = "Detected v2.txt")
#run_bands(mosaics, bands, 'plot', detected = "Detected v2.txt", fit = fit_func)
#run_bands(beams,'plot', detected = "Detected v2.txt", fit = fit_func)
#run_bands(['Pointings_CSV'],'plot', detected = "Detected v2.txt", fit = fit_func)
