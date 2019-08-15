
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
#from kinetics import load_dataframe


def func(S, Km, Vmax):
    V0 = (Vmax * S) / (Km + S)
    return V0

'''
def setup_data():
    # Probably don't need this function
    args = docopt.docopt(__doc__)
    if os.path.isfile(args['<bio96_metadata>']):
        data = bio96.load(args['<bio96_metadata>'],load_dataframe,{'well':'variable'})[0]
    elif os.path.isdir(args['<bio96_metadata>']):
        dataframes = []
        for f in os.listdir(args['<bio96_metadata>']):
            if f.endswith('.data'):
                subdata = bio96.load(os.path.join(args['<bio96_metadata>'],f),load_dataframe,{'well':'variable'})[0]
                #print(subdata)
                dataframes.append(subdata)
        data = pd.concat(dataframes,ignore_index=True)

    data = data.dropna()

    data['product'] = data['units'] * data['value']/data['conversion_factor']
    if 'enzyme_conc' in data:
        data['persecond'] = data['product'] / data['enzyme_conc']
    else:
        data['persecond'] = data['product']
    data['clicked_linear'] = False
    data['clicked_kinetics'] = False
    data['slope'] = 1
    return data
'''
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
ymax = max(data['slope']) + 0.1 * max(data['slope'])
ymin = min(data['slope']) - 0.1 * max(data['slope'])

# Colors
point_colors = cycle(['#C5523C', '#4BA0C1'])
line_colors = cycle(['#F06449', '#5BC3EB'])
area = 4*np.pi*3

def plot_stuff(graphWidth, graphHeight):
    # Plot
    f = plt.figure(figsize=(graphWidth/100.0,\
        graphHeight/100.0),dpi=100)
    axes = f.add_subplot(111)
    
    # raw data as scatter plot
    for group in parsed_data:
        axes.scatter(group[0], group[1], s=area, c=next(point_colors), alpha=0.5)
    for fit in fits:
        axes.plot(fit[0],fit[1], c=next(line_colors))

    axes.set_ylim([ymin, ymax])
    plt.show()

plot_stuff(800,600)
