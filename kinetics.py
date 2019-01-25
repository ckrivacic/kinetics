
"""

Usage:
    kinetics.py [options] <bio96_metadata>

Options (needs updating; some/most are depreciated and should instead be set
in the bio96 file):
    --conversion-factor NUMBER, -e NUMBER   [default: 1.0]
        Your units in the csv file will be divided by this.
        It is assumed that this will result in product concentration 
        in molar (M) units. Typically this would be an extinction
        coefficient with units M-1cm-1.

    --units INTEGER, -u INTEGER   [default: 0] 
        How many orders of magnitude to multiply the data by (after
        applying the conversion factor)

    --enzyme NUMBER, -n NUMBER   [default: 1.0]
        Enzyme concentration in the same units as substrate. 

    --start-row INTEGER, -s INTEGER   [default: 1]
        First row of header of data. 

    --end-row INTEGER, -d INTEGER
        Last row of data you want included. 

    --skip-rows LIST-LIKE, -k LIST-LIKE
        List or range of rows to skip. 
        Ex: 3 5 7:12
"""

import pandas as pd
from klab import docopt
from scipy import stats
import numpy as np
import scipy.optimize as opt
import os, bio96



def load_dataframe(path):

    args = docopt.docopt(__doc__)
    conversion_factor = float(args['--conversion-factor'])
    units = int(args['--units'])
    enzyme_conc = args['--enzyme']
    if enzyme_conc:
        enzyme_conc = float(enzyme_conc)
    else:
        enzyme_conc = 1.0

    start = int(args['--start-row']) - 1
    end = args['--end-row']

    skip = args['--skip-rows']

    """
    Prep the dataframe
    """

    if end:
        end = int(args['--end-row'])
        nrows = end - start - 1
        data = pd.read_csv(path,skiprows=start,nrows=nrows)

    else:
        data = pd.read_csv(path,skiprows=start)

    data = data.dropna(axis=1,how='all')
    #data = data.dropna()
    # Get time in units of seconds


    data['Time'] = data['Time'].str.split(':').apply(lambda x: int(x[0]) * 3600 + int(x[1])\
            * 60 + int(x[2]))
    for col in data:
        data[col] = pd.to_numeric(data[col])

    skiplist = []
    if skip:
        skip = skip.split(' ')
        for item in skip:
            if '-' in item:
                user_range = item.split('-')
                skiplist.extend(range(int(user_range[0])-start-2,int(user_range[1])
                    - start-1))
            else:
                skiplist.append(int(item)-start-2)

    data = data.drop(data.index[skiplist])
    #data = pd.melt(data, value_vars = 'value')


    for column in data:
        if column != 'Time':
            data[column] = data[column].multiply((10**units)/(conversion_factor * enzyme_conc))
    
    return pd.melt(data,id_vars=['Time'])


args = docopt.docopt(__doc__)

data = bio96.load(args['<bio96_metadata>'],load_dataframe,{'well':'variable'})[0]
data = data.dropna()

data['product'] = data['units'] * data['value']/data['conversion_factor']
data['persecond'] = data['product'] / data['enzyme_conc']
data['clicked_linear'] = False
data['clicked_kinetics'] = False
data['slope'] = 1


# The function we are trying to fit the data to
def func(S, Km, Vmax):
    V0 = (Vmax * S) / (Km + S)
    return V0

# Create class to store slope and concentration data
class regression_data(object):
    def __init__(self):
        self.data = {}
        self.slopes=[]
        self.concs=[]



import dash, dash_table
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from itertools import cycle
from flask_caching import Cache
from uuid import uuid4


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

CACHE_CONFIG= {
        'CACHE_TYPE': 'simple',
        }
cache= Cache()
cache.init_app(app.server, config=CACHE_CONFIG)


def serve_layout():
    max_time = data['Time'].max()
    session_id = str(uuid4())
    return html.Div([
        html.Div([

            dcc.Graph(id='linear-graph',
            ),

            dcc.RangeSlider(
                id='time-slider',
                min=data['Time'].min(),
                max=data['Time'].max(),
                value=[data['Time'].min(),data['Time'].max()],
                marks={(0):{'label':'0 min'},
                    str(data['Time'].max()): {'label':'{} min'.format(max_time)}},
                #marks={time: str(time) for time in data['Time']}),
                included=True,
                ),

            dcc.Graph(id='kinetics-graph'),

            dash_table.DataTable(id='kinetics-table',
                columns=(
                    [{'id':'enzyme','name':'Enzyme'},
                        {'id':'km','name':'Km'},
                        {'id':'kcat','name':'kcat'}]),
            )
            
            
        ]),

        html.Div(id='signal', style={'display':'none'}),
        html.Div(session_id,id='session_id', style={'display':'none'})

    ])

app.layout = serve_layout

colors = cycle(['rgb(57,106,177)','rgb(114,147,203)',
        'rgb(218,124,48)','rgb(225,151,76)',
        'rgb(62,150,81)','rgb(132,186,91)',
        'rgb(204,37,41)','rgb(211,94,96)',
        'rgb(83,81,84)','rgb(128,133,133)',
        'rgb(107,76,154)','rgb(144,103,167)',
        'rgb(146,36,40)','rgb(171,104,87)',
        'rgb(148,139,61)','rgb(204,194,16)'])

# Dictionary of dataframes indexed by session id. 
all_data = {}
kinetics_data = {}

"""
Lots of expensive calculations in the function below.
"""
@cache.memoize()
def update_kinetics_graph_global(value,session_id,clickData_kinetics=None):

    if session_id in kinetics_data:
        local_data = kinetics_data[session_id]
    else:
        return


    nonlinear_traces = []
    optimizedParameters_dict = {}
    # If we're here because a point was clicked on the kinetics graph,
    # update the local_data. 
    if clickData_kinetics:
        clickdat_indices = (local_data.index[(local_data['conc_uM']==clickData_kinetics['points'][0]['x']) &\
                (local_data['slope']==clickData_kinetics['points'][0]['y'])])
        for i in clickdat_indices:
            local_data['clicked_kinetics'].iloc[i] = not local_data['clicked_kinetics'][i]


    for name, group in local_data.groupby(['enzyme','clicked_kinetics']):
        current_color = next(colors)
        if name[1] == False:
            points = go.Scatter(
                    x = group['conc_uM'],
                    y = group['slope'],
                    mode='markers',
                    marker={'size':10,'color':current_color},
                    name=name[0])

            optimizedParameters,pcov = opt.curve_fit(func, group['conc_uM'],
                    group['slope'])
            optimizedParameters_dict[name[0]] = optimizedParameters
            xnew = np.linspace(min(group['conc_uM']),max(group['conc_uM']),100)
            ynew = func(xnew, optimizedParameters[0],optimizedParameters[1])

            line = go.Scatter(
                    x = xnew,
                    y = ynew,
                    marker={'color':next(colors)},
                    mode='lines',
                    name=name[0] + ' line')
            nonlinear_traces.append(points)
            nonlinear_traces.append(line)
        elif name[1] == True:
            points = go.Scatter(
                    x = group['conc_uM'],
                    y = group['slope'],
                    mode='markers',
                    marker={'symbol':'cross','size':10,'color':current_color},
                    name='Excluded points for ' + name[0])
            nonlinear_traces.append(points)
    return nonlinear_traces, optimizedParameters_dict


@cache.memoize()
def update_linear_graph_global(value,session_id,clickData_linear=None):
   
    if session_id in all_data:
        local_data = all_data[session_id]
    else:
        local_data = data.copy()
        all_data[session_id] = local_data
    groups = local_data.groupby(['enzyme', 'conc_uM', 'replicate'])
    traces = []
    x_min=value[0]
    x_max=value[1]

    # Keep track of slopes so we can make a smaller dataframe which is
    # quicker to fit to. 
    rows_to_add_to_kinetics_df = []

    for name,group in groups:
        enzyme = name[0]
        conc = name[1]
        replicate = name[2]

        df = group[(group['Time'] >= x_min) & (group['Time'] <= x_max)]

        xi = df['Time']
        yi = df['persecond']
        trace1 = go.Scatter(
            x = xi,
            y = yi,
            mode='markers',
            marker={'size':7,'color':next(colors)},
            name=str(enzyme) + ', ' + str(conc) + ' uM '+ str(replicate)
        )
        traces.append(trace1)

        slope, intercept, r_value, p_value, std_err = stats.linregress(xi,yi)
        line = slope * xi + intercept

        rows_to_add_to_kinetics_df.append([enzyme, conc, slope,
                replicate])
        
        trace2 = go.Scatter(
            x = xi,
            y = line,
            mode='lines',
            marker={'color':next(colors)},
            name=str(enzyme) + ', ' + str(conc) + ' uM '+ str(replicate) + 'line')
        traces.append(trace2)

        local_data.loc[(local_data['enzyme']==enzyme) & \
        (local_data['conc_uM']==conc) & (local_data['replicate']==replicate),'slope'] = slope

    dict_list = []
    for row in rows_to_add_to_kinetics_df:
        dict1 = {}
        dict1['enzyme'] = row[0]
        dict1['conc_uM'] = row[1]
        dict1['slope'] = row[2]
        dict1['replicate'] = row[3]

        dict_list.append(dict1)
    kinetics_df = pd.DataFrame(dict_list)
    kinetics_df['clicked_kinetics'] = False

    kinetics_data[session_id] = kinetics_df


    update_kinetics_graph_global(value,session_id)
    return traces



"""
Signal callback
"""
@app.callback(
    dash.dependencies.Output('signal','children'),
    [dash.dependencies.Input('time-slider','value'),
        ])
def update_graphs(value):
    #global_store(value,session_id)
    return value



"""
Update linear graph
"""

@app.callback(
        dash.dependencies.Output('linear-graph','figure'),
        [dash.dependencies.Input('signal','children'),
            dash.dependencies.Input('session_id','children')])
def update_linear_graph(value,session_id):
    traces = update_linear_graph_global(value,session_id)
    return{
        'data': traces,
        'layout': go.Layout(
            xaxis={
                'title':'Time (s)',
                'type': 'linear'
            },
            yaxis={
                'title':'product',
                'type': 'linear'
            },
            margin={'l':40,'b':40,'t':100,'r':0},
            hovermode='closest',
            title='Time courses'
        )
    }


"""
Update kinetics graph
"""

@app.callback(
    dash.dependencies.Output('kinetics-graph','figure'),
    [dash.dependencies.Input('signal','children'),
        dash.dependencies.Input('session_id','children'),
        dash.dependencies.Input('kinetics-graph','clickData')])
def update_kinetics_graph(value,session_id,clickData):
    print(clickData)
    nonlinear_traces = update_kinetics_graph_global(value,session_id,clickData_kinetics=clickData)[0]
    optimizedParameters = update_kinetics_graph_global(value,session_id,clickData_kinetics=clickData)[1]
    
    return{
        'data':nonlinear_traces,
        'layout':go.Layout(
            xaxis={
                'title':'Concentration (uM)'
            },
            yaxis={'title':'V0'},
            margin={'l':40,'b':40,'t':100,'r':10},
            hovermode='closest',
            title='Michaelis-Menten Kinetics'
         )
    }

"""
Update kinetics table
"""

@app.callback(
    dash.dependencies.Output('kinetics-table','data'),
    [dash.dependencies.Input('signal','children'),
        dash.dependencies.Input('session_id','children')])
def update_data_table(value,session_id):
    optimizedParameters = update_kinetics_graph_global(value,session_id)[1]
    table_data = []
    for enzyme in optimizedParameters:
        table_data.append({'enzyme':enzyme, 'km':optimizedParameters[enzyme][0]\
                , 'kcat':optimizedParameters[enzyme][1]})
    return table_data


if __name__ == '__main__':
    app.run_server(debug=True)