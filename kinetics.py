
"""

Usage:
    kinetics.py [options] (--load <pklfile> | <bio96_metadata>...)

Options (needs updating; some/most are depreciated and should instead be set
in the bio96 file):
    --conversion-factor NUMBER, -e NUMBER   [default: 1.0]
        Your units in the csv file will be divided by this.
        It is assumed that this will result in product concentration 
        in molar (M) units. Typically this would be an extinction
        coefficient with units M-1cm-1.

    --correction NUMBER, -c NUMBER [default:1.0]
        Take product row and divide it by this number. This is for going from cuvettes to wells. 

    --raw
        Plot raw data only. 

    --zero, -z
        Subtract zero substrate

    --zerokinetics, -i

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
import scipy, plotly
from scipy import stats
import numpy as np
import pickle as pkl
import scipy.optimize as opt
import os, bio96, sys

args = docopt.docopt(__doc__)

def convert_time(x):
    seconds = float(x[0]) * 3600 + float(x[1]) * 60 + float(x[2])
    return seconds

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
    
    data['Time'] = data['Time'].str.split(':').apply(lambda x:\
            convert_time(x))
            #float(x[0]) * 3600 + float(x[1])\
            #* 60 + float(x[2]))

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

def get_data():
    if os.path.isdir(args['<bio96_metadata>'][0]):
        dataframes = []
        for f in os.listdir(args['<bio96_metadata>']):
            if f.endswith('.data'):
                subdata = bio96.load(os.path.join(args['<bio96_metadata>'],f),load_dataframe,{'well':'variable'})
                #print(subdata)
                dataframes.append(subdata)
        data = pd.concat(dataframes,ignore_index=True)
    else:
        dataframes = []
        for f in args['<bio96_metadata>']:
            subdata = bio96.load(f, load_dataframe,{'well':'variable'})
            dataframes.append(subdata)
        data = pd.concat(dataframes, ignore_index=True)

    data = data.dropna()
    raw = args['--raw']

    # Subtract the no-substrate well
    if args['--zero']:
        import time
        start_time = time.time()
        zero_df = data[data['conc_uM']==0]
        zgroups = zero_df.groupby(['enzyme','date','replicate'])
        groupdict = {}
        for name,group in zgroups:
            if name[0] not in groupdict:
                groupdict[name[0]]={}
            if name[1] not in groupdict[name[0]]:
                groupdict[name[0]][name[1]]={}
            groupdict[name[0]][name[1]][name[2]]=group
        for index, row in data.iterrows():
            enzyme = row['enzyme']
            date = row['date']
            replicate = row['replicate']
            #print(replicate)
            conc = row['conc_uM']
            timept = row['Time']
            zdf = groupdict[enzyme][date][replicate]
            zero = zdf.loc[zdf['Time']==timept]['value'].tolist()[0]
            """
            zero = zero_df.loc[(zero_df['enzyme']==enzyme) & \
                    (zero_df['date']==date) & \
                    (zero_df['replicate']==replicate) & \
                    (zero_df['conc_uM']==conc) & \
                    (zero_df['Time']==time),
                    ['conc_uM']]
            """
            row['value'] = row['value'] - zero
            data.iloc[index] = row
    #print('time to load:', time.time() - start_time)


    if raw:
        data['product'] = data['value']
    else:
        data['product'] = data['units'] * data['value']/data['conversion_factor']

    if args['--correction']:
        correction = float(args['--correction'])
        data['product'] = data['product'] * correction


    if 'enzyme_conc' in data and not raw:
        data['persecond'] = data['product'] / data['enzyme_conc']
    else:
        data['persecond'] = data['product']
    if 'enzyme_conc' not in data:
        data['enzyme_conc'] = 1
    data['clicked_linear'] = False
    data['clicked_kinetics'] = False
    data['slope'] = 1
    data['min_time'] = min(data['Time'])
    data['max_time'] = max(data['Time'])
    data['shown'] = True

    return data

def load_session(pklfname):
    with open(pklfname, 'rb') as pklfile:
        analysis_data = pkl.load(pklfile)
    return analysis_data

kinetics = None
if args['<pklfile>']:
    analysis_data = load_session(args['<pklfile>'])
    data = analysis_data['local_data']
    kinetics = analysis_data['kinetics_df']
elif args['<bio96_metadata>']:
    data = get_data()
else:
    print('No data given. Either provide a session via --load or point to a bio96 metadata file.')

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
from itertools import cycle,islice
from flask_caching import Cache
from uuid import uuid4


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

CACHE_CONFIG= {
        'CACHE_TYPE': 'simple',
        }
cache= Cache()
cache.init_app(app.server, config=CACHE_CONFIG)


# Dictionary of dataframes indexed by session id. 
all_data = {}
kinetics_data = {}


def serve_layout(data):
    max_time = data['Time'].max()
    session_id = str(uuid4())

    if kinetics is not None:
        kinetics_data[session_id] = kinetics
        all_data[session_id] = data

    enzyme_list = []
    for enzyme in set(data[data['shown']==True]['enzyme']):
        enzyme_list.append({'label':enzyme,'value':enzyme})

    conc_list = []
    for conc in set(data[data['shown']==True]['conc_uM']):
        conc_list.append({'label':str(conc) + ' uM','value':conc})

    replicate_list = []
    for rep in set(data[data['shown']==True]['replicate']):
        replicate_list.append({'label':'Replicate ' + str(rep), 'value':rep})

    date_list = []
    for date in set(data[data['shown']==True]['date']):
        date_list.append({'label':date,'value':date})

    return html.Div([
        html.Div([

            dcc.Dropdown(
                id='enzyme-dropdown',
                options=enzyme_list,
                placeholder='Filter by enzyme...',
                multi=True
            ),

            dcc.Dropdown(
                id='conc-dropdown',
                options=conc_list,
                placeholder='Filter by concentration...',
                multi=True
            ),

            dcc.Dropdown(
                id='rep-dropdown',
                options=replicate_list,
                placeholder='Filter by replicate...',
                multi=True
            ),

            dcc.Dropdown(
                id='date-dropdown',
                options=date_list,
                placeholder='Filter by date...',
                multi=True
            ),

            html.Button(id='filter-button',n_clicks=0,children='Filter'),

            dcc.Graph(id='linear-graph',
            ),
            dcc.RadioItems(
                id='fit-type',
                options=[
                    {'label':'Linear','value':'linear'},
                    {'label':'Exponential','value':'exponential'}
                ],
                value='linear'
                ),
            dcc.Input(id='min-time',type='number',value=data['Time'].min()),
            dcc.Input(id='max-time',type='number',value=data['Time'].max()),
            dcc.Input(id='max-percent-substrate',type='number',value=100),

            html.Button(id='submit-button',n_clicks=0,children='Submit'),

            html.Button(id='save-button',n_clicks=0,children='Save analysis progress'),

            dcc.Input(id='p0-c',type='number',value='50'),
            dcc.Input(id='p0-s0',type='number',value='30'),
            dcc.Input(id='p0-k',type='number',value='0.01'),

            dcc.Graph(id='kinetics-graph'),

            dash_table.DataTable(id='kinetics-table',
                columns=(
                    [{'id':'enzyme','name':'Enzyme'},
                        {'id':'km','name':'Km'},
                        {'id':'km stdev','name':'Km std dev errors'},
                        {'id':'kcat','name':'kcat'},
                        {'id':'kcat stdev','name':'kcat std dev errors'}]),
            ),

            dcc.Textarea(
                    id='log',
                    style={'width':'100%'},
                    disabled='True',
                    value='Log info goes here'
                )
            
            
            ],style={'padding':50}),

        html.Div(session_id,id='session_id',
            style={'display':'none','padding':10})

        ],style={'padding':100})

app.layout = serve_layout(data)

colors = ['rgb(57,106,177)','rgb(114,147,203)',
        'rgb(218,124,48)','rgb(225,151,76)',
        'rgb(62,150,81)','rgb(132,186,91)',
        'rgb(204,37,41)','rgb(211,94,96)',
        'rgb(83,81,84)','rgb(128,133,133)',
        'rgb(107,76,154)','rgb(144,103,167)',
        'rgb(146,36,40)','rgb(171,104,87)',
        'rgb(148,139,61)','rgb(204,194,16)']



"""
Lots of expensive calculations in the function below.
"""
#@cache.memoize()
def update_kinetics_graph_global(value,session_id,clickData_kinetics=None):

    if session_id in kinetics_data:
        local_data = kinetics_data[session_id]
    else:
        return


    nonlinear_traces = []
    optimizedParameters_dict = {}
    err_dict = {}
    # If we're here because a point was clicked on the kinetics graph,
    # update the local_data. 
    if clickData_kinetics:
        clickdat_indices = (local_data.index[(local_data['conc_uM']==clickData_kinetics['points'][0]['x']) &\
                (local_data['slope']==clickData_kinetics['points'][0]['y'])])
        for i in clickdat_indices:
            local_data.at[i,'clicked_kinetics'] = not local_data['clicked_kinetics'][i]
    kinetics_data[session_id] = local_data
    color_cycle = cycle(colors)
    current_color = colors[0]

    if args['--zerokinetics']:
        try:
            zero_df = local_data[local_data['conc_uM']==0]
            zgroups = zero_df.groupby(['enzyme','date'])
            groupdict = {}
            for name,group in zgroups:
                if name[0] not in groupdict:
                    groupdict[name[0]]={}
                if name[1] not in groupdict[name[0]]:
                    groupdict[name[0]][name[1]]={}
                groupdict[name[0]][name[1]]=group
            for index, row in local_data.iterrows():
                enzyme = row['enzyme']
                date = row['date']
                replicate = row['replicate']
                conc = row['conc_uM']
                slope = row['slope']

                try:
                    zdf = groupdict[enzyme][date]
                    zero = zdf[zdf['clicked_kinetics']==False]['slope'].tolist()
                    #zero= zdf['slope'].tolist()
                    zavg = sum(zero) / len(zero)

                    row['slope'] = row['slope'] - zavg
                    local_data.iloc[index] = row
                except:
                    break
        except:
            print("skipping zeroing")
    for name, group in local_data.groupby(['date','enzyme','clicked_kinetics'],sort=True):
        if name[2] == False:
            current_color = next(color_cycle)
            points = go.Scatter(
                    x = group['conc_uM'],
                    y = group['slope'],
                    mode='markers',
                    marker={'size':10,'color':current_color},
                    name=name[1] + '_' + name[0])
            optimizedParameters,pcov = opt.curve_fit(func, group['conc_uM'], group['slope'])
            err_dict[name[1]+name[0]] = np.sqrt(np.diag(pcov))
            optimizedParameters_dict[name[1]+name[0]] = optimizedParameters
            xnew = np.linspace(min(group['conc_uM']),max(group['conc_uM']),100)
            ynew = func(xnew, optimizedParameters[0],optimizedParameters[1])

            line = go.Scatter(
                    x = xnew,
                    y = ynew,
                    marker={'color':next(color_cycle)},
                    mode='lines',
                    name=name[1] + '_' + name[0] + ' line')
            nonlinear_traces.append(points)
            nonlinear_traces.append(line)
        elif name[2] == True:
            points = go.Scatter(
                    x = group['conc_uM'],
                    y = group['slope'],
                    mode='markers',
                    marker={'symbol':'cross','size':10,'color':current_color},
                    name='Excluded points for ' + name[1] + '_' + name[0])
            nonlinear_traces.append(points)
    return nonlinear_traces, optimizedParameters_dict, err_dict


#@cache.memoize()
def update_linear_graph_global(value,session_id,fit_type='linear',clickData_linear=None,max_percent_substrate=None,
        filters={'enzyme':[],'conc_uM':[],'replicate':[],'date':[]},update_time=True,update_linear=True,p0_c=50,p0_s0=30,p0_k=0.01):

    if session_id in all_data:
        local_data = all_data[session_id]
    else:
        local_data = data.copy()
        all_data[session_id] = local_data

    if filters:
        for f in filters:
            if not filters[f]: # is the list empty?
                filters[f] = list(set(local_data[f]))

    local_data['shown'] = False
    local_data.loc[(local_data['enzyme'].isin(filters['enzyme'])) & \
            (local_data['date'].isin(filters['date'])) & \
            (local_data['replicate'].isin(filters['replicate'])) & \
            (local_data['conc_uM'].isin(filters['conc_uM'])),'shown'] = True
    if update_time:
        local_data.loc[local_data['shown'] == True,'min_time'] = value[0]
        local_data.loc[local_data['shown'] == True,'max_time'] = value[1]
    if clickData_linear:
        clickdat_indices = (local_data.index[(local_data['Time']==clickData_linear['points'][0]['x']) &\
                (local_data['product']==clickData_linear['points'][0]['y'])])
        for i in clickdat_indices:
            local_data.at[i,'clicked_linear'] = not local_data['clicked_linear'][i]
    groups = local_data.groupby(['date','enzyme','conc_uM', 'replicate','clicked_linear','shown'],sort=True)
    traces = []

    # Keep track of slopes so we can make a smaller dataframe which is
    # quicker to fit to. 
    rows_to_add_to_kinetics_df = []

    colors_cycle = cycle(colors)
    current_color = colors[0]
    for name,group in groups:
        date = name[0]
        enzyme = name[1]
        conc = name[2]
        replicate = name[3]
        clicked = name[4]
        shown = name[5]

        if max_percent_substrate and conc != 0:
            min_product = min(group['product'])
            group = group[(group['product'] - min_product) / group['conc_uM'] <= max_percent_substrate/100]
        df = group[(group['Time'] >= group['min_time']) & (group['Time'] <= group['max_time'])]

        if clicked == False and shown == True:
            current_color = next(colors_cycle)
            xi = df['Time']
            yi = df['product']
            trace1 = go.Scatter(
                x = xi,
                y = yi,
                mode='markers',
                marker={'size':7,'color':current_color},
                name=str(enzyme) + ', ' + str(conc) + ' uM '+ str(replicate) + '_' + str(date)
            )
            traces.append(trace1)

            slope = None
            line = None
            skip = False
            plateau = None

            enzyme_conc = list(set(group['enzyme_conc']))[0]
            if fit_type == 'linear': #or conc == 0 or enzyme == 'None':
                try:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(xi,yi)
                    line = slope * xi + intercept
                    if len(xi) > 3:
                        rows_to_add_to_kinetics_df.append([enzyme, conc, slope/enzyme_conc,
                                replicate,date])
                    if len(xi)==1:
                        skip = True
                except:
                    print('No fit found for ' + enzyme + ' replicate ' + str(replicate) + ' conc. ' + str(conc))
                    skip = True
            elif fit_type == 'exponential':
                #if conc != 0:
                try:
                
                #cguess = conc/(10*450)
                #cguess = 80
                    print('Fitting using the following initial parameters for A0, c, and k:')
                    print(p0_s0, p0_c, p0_k)

                    s0,c, slope = scipy.optimize.curve_fit(lambda t, s0,c, k: c - s0 * np.exp(-k * t), xi, yi,p0=(p0_s0,p0_c,p0_k),maxfev=2000)[0]
                    line = (c - s0 * np.exp(-slope * xi))
                    print('Fit ' + enzyme + ' replicate ' + str(replicate) + ' conc. ' + str(conc) + ' with the following parameters:')
                    print('A0: ' + str(s0))
                    print('c: ' + str(c))
                    print('k: ' + str(slope))
                    plateau = c
                    v0 = (slope * s0) / enzyme_conc
                    rows_to_add_to_kinetics_df.append([enzyme, conc, v0,
                            replicate,date])
                except:
                    print('No fit found for ' + enzyme + ' replicate ' + str(replicate) + ' conc. ' + str(conc))
                    skip = True

            if not skip: 
                trace2 = go.Scatter(
                    x = xi,
                    y = line,
                    mode='lines',
                    marker={'color':next(colors_cycle)},
                    name=str(enzyme) + ', ' + str(conc) + ' uM '+ str(replicate)  + '_' + str(date) + 'line')
                traces.append(trace2)

                if fit_type == 'linear':
                    local_data.loc[(local_data['enzyme']==enzyme) & \
                    (local_data['conc_uM']==conc) & (local_data['replicate']==replicate) & \
                    (local_data['date']==date),'slope'] = slope
                elif fit_type == 'exponential':
                    local_data.loc[(local_data['enzyme']==enzyme) & \
                    (local_data['conc_uM']==conc) & (local_data['replicate']==replicate) & \
                    (local_data['date']==date),'slope'] = plateau
            else:
                local_data.loc[(local_data['enzyme']==enzyme) & \
                (local_data['conc_uM']==conc) & (local_data['replicate']==replicate) & \
                (local_data['date'] == date),'slope'] = 0

        elif clicked == True and shown == True:
            xi = df['Time']
            yi = df['product']
            points = go.Scatter(
                    x = xi,
                    y = yi,
                    mode='markers',
                    marker={'symbol':'cross','size':10,'color':current_color},
                    name='Excluded points for ' + str(enzyme) + ', ' + str(conc) + ' uM '+ str(replicate)  + '_' + str(date)
                    )
            traces.append(points)
            
    if update_linear == True:
        dict_list = []
        for row in rows_to_add_to_kinetics_df:
            dict1 = {}
            dict1['enzyme'] = row[0]
            dict1['conc_uM'] = row[1]
            dict1['slope'] = row[2]
            dict1['replicate'] = row[3]
            dict1['date'] = row[4]

            dict_list.append(dict1)
        kinetics_df = pd.DataFrame(dict_list)

        kinetics_df['clicked_kinetics'] = False

        kinetics_data[session_id] = kinetics_df


    #update_kinetics_graph_global(value,session_id)
    return traces


general_data = {}

@app.callback([dash.dependencies.Output('log','value')],
        [dash.dependencies.Input('save-button','n_clicks')],
        [dash.dependencies.State('session_id','children'),
            dash.dependencies.State('log','value')]
        )
def save_progress(n_clicks,session_id,logtext):
    analysis_data = {}
    analysis_data['local_data'] = all_data[session_id]
    analysis_data['kinetics_df'] = kinetics_data[session_id]
    infile = args['<bio96_metadata>']
    if infile:
        if os.path.isfile(infile[0]):
            outfname = infile[0] + str(session_id) + '.pkl'
        else:
            outfname = os.path.join(infile,session_id + '.pkl')
    else:
        outfname = args['<pklfile>']

    newtext = 'Saving to ' + outfname
    
    with open(outfname,'wb') as outfile:
        pkl.dump(analysis_data,outfile)
    return([logtext + '\n' + newtext])





"""
Update linear graph
"""


@app.callback(
        [dash.dependencies.Output('linear-graph','figure'),
            dash.dependencies.Output('kinetics-graph','figure'),
            dash.dependencies.Output('kinetics-table','data'),
            ],        
        [dash.dependencies.Input('submit-button','n_clicks'),
            dash.dependencies.Input('filter-button','n_clicks'),
            dash.dependencies.Input('linear-graph','clickData'),
            dash.dependencies.Input('kinetics-graph','clickData'),
        ],
        [
            dash.dependencies.State('session_id','children'),
            dash.dependencies.State('max-percent-substrate','value'),
            dash.dependencies.State('min-time','value'),
            dash.dependencies.State('max-time','value'),
            dash.dependencies.State('fit-type','value'),
            dash.dependencies.State('enzyme-dropdown','value'),
            dash.dependencies.State('conc-dropdown','value'),
            dash.dependencies.State('rep-dropdown','value'),
            dash.dependencies.State('date-dropdown','value'),
            dash.dependencies.State('p0-c','value'),
            dash.dependencies.State('p0-s0','value'),
            dash.dependencies.State('p0-k','value')
            ])
def update_linear_graph(n_clicks_submit,n_clicks_filter,clickData_linear,clickData_kinetics,session_id,percent_substrate_value,min_time,max_time,fit_type,
        enzyme_filter,conc_filter,rep_filter,date_filter,p0_c,p0_s0,p0_k):
    """
    Update linear graph
    """
    timeslider_value = [min_time,max_time]
    if not session_id in general_data:
        general_data[session_id] = {'min-time':min_time,'max-time':max_time,'percent_substrate_value':percent_substrate_value,'clickData_linear':clickData_linear,'clickData_kinetics':clickData_kinetics,'n_clicks_submit':n_clicks_submit,'n_clicks_filter':n_clicks_filter,'p0_c':p0_c,'p0_s0':p0_s0,'p0_k':p0_k}

    update_time = False 
    if n_clicks_submit==general_data[session_id]['n_clicks_submit']:
        update_time = False
    else:
        update_time = True
   
    update_linear = True
    if clickData_kinetics != general_data[session_id]['clickData_kinetics']:
        update_linear = False
        general_data[session_id]['clickData_kinetics'] = clickData_kinetics
    linear_traces = update_linear_graph_global(timeslider_value,session_id,fit_type,clickData_linear = clickData_linear, max_percent_substrate=percent_substrate_value,
            filters={'enzyme':enzyme_filter,'conc_uM':conc_filter,'replicate':rep_filter,'date':date_filter},update_time=update_time,update_linear=update_linear,p0_c=float(p0_c),\
                    p0_s0=float(p0_s0),p0_k=float(p0_k))
    linear_output = {
        'data': linear_traces,
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
    Update nonlinear graph
    """
    try:
        kinetics_update = update_kinetics_graph_global(timeslider_value,session_id,clickData_kinetics=clickData_kinetics)
        nonlinear_traces = kinetics_update[0]
        optimizedParameters = kinetics_update[1]
        errorParameters = kinetics_update[2]
    except:
        nonlinear_traces = None
        optimizedParameters = None
        errorParameters = None

    nonlinear_output = {
        'data':nonlinear_traces,
        'layout':go.Layout(
            xaxis={
                'title':'Concentration (uM)'
            },
            yaxis={'title':'','position':0},
            margin={'l':40,'b':40,'t':100,'r':10},
            hovermode='closest',
            title='Michaelis-Menten Kinetics'
         )
    }

    """
    Update table
    """
    table_data = []
    if optimizedParameters and errorParameters:
        for enzyme in optimizedParameters:
            table_data.append({'enzyme':enzyme, 'km':optimizedParameters[enzyme][0]\
                    ,'km stdev':errorParameters[enzyme][0], 'kcat':optimizedParameters[enzyme][1],\
                    'kcat stdev':errorParameters[enzyme][1]})

    general_data[session_id] = {'min-time':min_time,'max-time':max_time,'percent_substrate_value':percent_substrate_value,'clickData_linear':clickData_linear,'clickData_kinetics':clickData_kinetics,'n_clicks_submit':n_clicks_submit,'n_clicks_filter':n_clicks_filter,'p0_c':p0_c,'p0_s0':p0_s0,'p0_k':p0_k}
    return linear_output, nonlinear_output, table_data



if __name__ == '__main__':

    app.run_server(debug=True,host='0.0.0.0',port=8080)
