
"""

Usage:
    kinetics.py [options] (<csv>...)

"""
import numpy as np
import pandas as pd
import scipy.optimize as opt
from klab import docopt
from itertools import cycle
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.pyplot as pyplot


def func(S, Km, Vmax):
    V0 = (Vmax * S) / (Km + S)
    return V0

def import_data():
    args = docopt.docopt(__doc__)
    dataframes = []
    for csv in args['<csv>']:
        subdata = pd.read_csv(csv)
        dataframes.append(subdata)
    return pd.concat(dataframes, ignore_index=True)


# Grab data
data = import_data()
# Get rid of excluded points
data = data[data['clicked_kinetics']==False]
# Create data
parsed_data = []
fits = []
for name, group in data.groupby(['date','enzyme']):
    x = group['conc_uM']
    y = group['slope']
    parsed_data.append((x,y))

    # Generate parameters
    optimizedParameters, pcov = opt.curve_fit(func, x, y)
    xModel = np.linspace(min(x),max(x),num=500)
    yModel = func(xModel, *optimizedParameters)
    fits.append((xModel, yModel))

# Axes limits
ymax = max(data['slope']) + 0.05 * max(data['slope'])
ymin = min(data['slope']) - 0.05 * max(data['slope'])
xmax = max(data['conc_uM']) + 0.05 * max(data['conc_uM'])
xmin = min(data['conc_uM']) - 0.05 * max(data['conc_uM'])

# Colors
point_colors = cycle(['#C5523C', '#4BA0C1'])
line_colors = cycle(['#F06449', '#5BC3EB'])
area = 4*np.pi*3

def plot_stuff(graphWidth, graphHeight):
    font = {'family' : 'normal',
            'size' : 12}
    plt.rc('font', **font)

    # Plot
    f = plt.figure(figsize=(graphWidth/100.0,\
        graphHeight/100.0),dpi=100)
    axes = f.add_subplot(111)
    
    # raw data as scatter plot
    for group in parsed_data:
        axes.scatter(group[0], group[1], s=area, c=next(point_colors), alpha=0.5)
    # plot the fits
    for fit in fits:
        axes.plot(fit[0],fit[1], c=next(line_colors))

    pyplot.hlines(0,xmin,xmax,linestyles='dashed',label='')

    # Label parameters
    axes.set_ylim([ymin, ymax])
    axes.set_xlim([xmin, xmax])
    axes.tick_params(width=2,length=8)
    axes.tick_params(axis='y',which='both',labelleft=True,labelright=False,right=True)
    axes.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    axes.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    axes.xaxis.set_tick_params(direction='in')
    axes.yaxis.set_tick_params(direction='in')

    plt.xlabel('[5(10) estrene-dione], μM')
    plt.ylabel('V$_{0}$, s$^{-1}$')

    plt.show()

plot_stuff(800,600)
