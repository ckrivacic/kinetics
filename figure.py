import numpy as np
from klab import docopt
import matplotlib.pyplot as plt
from kinetics import load_dataframe

def func(S, Km, Vmax):
    V0 = (Vmax * S) / (Km + S)
    return V0

def setup_data:

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

# Create data
N = 500
x = np.random.rand(N)
y = np.random.rand(N)
colors = (0,0,0)
area = np.pi*3

def plot_stuff(graphWidth, graphHeight):
    # Plot
    f = plt.figure(figsize=(graphWidth/100.0,\
        graphHeight/100.0),dpi=100)
    axes = f.add_subplot(111)
    
    # raw data as scatter plot
    axes.scatter(x, y, s=area, c=colors, alpha=0.5)
    plt.title('Scatter plot pythonspot.com')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

plot_stuff(800,600)
