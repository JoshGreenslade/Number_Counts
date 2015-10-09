# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 12:40:22 2015

@author: JG1714

This is a section of code that will estimate the differential and integral 
number counts of sources from a catalogue.

In general the catalogue should only need the list of fluxes.

Before passing the list, 
"""

# =============================================================================
# Integral counts of the blank field number sources (schecter function)

nq = 1599
sp = 3.3
alpha = -2.0

# Define the schecter function
schec = lambda s: (nq / sp) * s * (np.power((s / sp),alpha)) * np.exp(-s / sp)

fluxes = np.array([])
n = 1
ma = 30
mi = 1

integral_counts = integrate.quad(schec, mi, ma)
s_max = schec(mi)
s_min = schec(ma)
while np.size(fluxes) < int(integral_counts[0]*0.2):
    n += 1
    test =  np.random.rand()*np.abs(s_max-s_min) + s_min
    test2 = np.random.rand()*np.abs(ma-mi)+mi
    
    
    
    if test <= schec(test2):
        fluxes = np.append(fluxes, test2)
        print(n)


fluxes = fluxes #### !!! Replace this line with the flux values !!! ###
area = 0.2 ### !!! Replace this line with the area (in square degrees) of survey.



####################### Integral Number Counts #############################

"""
This section will calculate the integral number counts from the catalogue.
"""

# Input parameters to the system

flux_minimum = 2   # The minimum flux to count from
flux_maximum = 20    # The maximum flux to count up to
n_bins = 30          # The number of bins to use



"""
For each bin, calculate the number of sources in it
"""

# Set up the bins   
bin_edges = np.linspace(flux_minimum, 
                        flux_maximum, 
                        n_bins + 1)  # As we will end up with n-1 bin centres
bin_centres = ( bin_edges + np.roll(bin_edges, -1) )/2
bin_centres = bin_centres[:-1]
bin_counts = np.zeros(n_bins) # n_bins - 1 as only values between values


# For each bin, add all the counts of sources within that bin

for index, Bin in enumerate(bin_edges[:-1]):
    
    n_sources = np.sum(fluxes >= Bin)
    bin_counts[index] = n_sources
    

# The errors on each bin are poisson in nature, so sqrt for error

bin_errors = np.sqrt(bin_counts)


# Convert these to the raw bin counts, so we can modify them later yet retain
# their values.

raw_Ibin_counts = bin_counts
raw_Ibin_errors = bin_errors

raw_int_counts = bin_counts/area
raw_int_errors = bin_errors/area
raw_int_centres = bin_centres

"""
We now have a simple count of the number of sources that have fluxes greater
than each value in bin_edges, and an estimate of the poisson error on each of 
them.

"""




####################### Differential Number Counts #############################

"""
This section will calculate the differential number counts from the catalogue.
Very similar to the previous section
"""


bin_widths = (np.roll(bin_edges, -1) - bin_edges)
bin_widths = bin_widths[:-1]

bin_counts = np.zeros(n_bins) # n_bins - 1 as only values between values


# For each bin, add all the counts of sources within that bin

for index, Bin in enumerate(bin_edges[:-1]):
    
    source_selection = (fluxes >= bin_edges[index]) * (fluxes < bin_edges[index+1]) 
    n_sources = np.sum(source_selection)
    bin_counts[index] = n_sources
    

# The errors on each bin are poisson in nature, so sqrt for error

bin_errors = np.sqrt(bin_counts)


# Convert these to the raw bin counts, so we can modify them later yet retain
# their values.

raw_Dbin_counts = bin_counts/bin_widths
raw_Dbin_errors = bin_errors/bin_widths

raw_dif_centres = bin_centres
raw_dif_counts = raw_Dbin_counts/area
raw_dif_errors = raw_Dbin_errors/area
