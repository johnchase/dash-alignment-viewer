"""Alignmet viewer application
"""
import os
import dash

import numpy as np

import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

from plotly import tools
from util import (alignment_layout, get_msa_order, get_dimensions,
                  parse_seq_object, parse_sequences)
from style import BASE_DIC, UPLOAD_BUTTON, COLOR_DIC


app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(children=[
    html.H1(children='Dash Alignment Viewer'),
    html.Div([

        html.Div([dcc.Upload(html.Button('Upload File',
                                         id='upload-button',
                                         n_clicks_timestamp='0',
                                         style=UPLOAD_BUTTON),
                             id='upload_data',
                             style={'display': 'inline-block'}),
                 ]),

    html.Div([

        dcc.Checklist(id='reference_layout',
                      options=[{'label': 'Reference Layout',
                                'value': True}],
                      values=[True],
                      style={'marginRight': 20,
                             'display': 'inline-block'}),


        dcc.RadioItems(id='layout-type',
                       options=[{'label': i, 'value': i} for i in ['Block',
                                                                   'Letter']],
                       value='Block',
                       style={'marginRight': 20,
                              'display': 'inline-block'}),

            dcc.Dropdown(id='parent-seq',
                style={'width': '200px',
                    'marginRight': 20,
                       'display': 'inline-block'}),

            dcc.Dropdown(id='color-palette',
                options=[{'label': key, 'value': key} for key, value 
                         in COLOR_DIC.items()],
                value='Dark2',
                style={'width': '200px', 'display': 'inline-block'}),
          
           html.Div([dcc.Dropdown(id='sample-data',
                                  options=[{'label': 'Example 1', 
                                            'value': os.path.join('data','example1.msa')},
                                           {'label': 'Example 2',
                                            'value': os.path.join('data',
                                                'msa10.fna')}],
                                  style={'width': '200px'})],
                id='sample-data-div',
                n_clicks_timestamp='0',
                style={'display': 'inline-block'})
        
        ]),

        dcc.Graph(
            id='alignment',
            config={
                'displayModeBar': False})])])


@app.callback(
    dash.dependencies.Output('parent-seq', 'options'),
    [dash.dependencies.Input('upload_data', 'contents'),
     dash.dependencies.Input('upload-button', 'n_clicks_timestamp'),
     dash.dependencies.Input('sample-data', 'value'),
     dash.dependencies.Input('sample-data-div', 'n_clicks_timestamp')])
def get_sequence_names(upload_object, upload_timestamp,
                       sample_object, sample_timestamp):
    
    
    if int(upload_timestamp) > int(sample_timestamp):
        seq_object = upload_object
        if seq_object is None:
            return ''
        seq_lines = parse_seq_object(seq_object)
    else:
        if sample_object is None:
            return ''
        seq_lines = sample_object
    names, _, _ = parse_sequences(seq_lines)
    return [{'label': label, 'value': label} for label in names]


@app.callback(
    dash.dependencies.Output('alignment', 'figure'),
    [dash.dependencies.Input('layout-type', 'value'),
     dash.dependencies.Input('parent-seq', 'value'),
     dash.dependencies.Input('upload_data', 'contents'),
     dash.dependencies.Input('upload-button', 'n_clicks_timestamp'),
     dash.dependencies.Input('sample-data', 'value'),
     dash.dependencies.Input('sample-data-div', 'n_clicks_timestamp'),
     dash.dependencies.Input('color-palette', 'value'),
     dash.dependencies.Input('reference_layout', 'values')])
def create_alignment(layout, reference_name,
                     upload_object, upload_timestamp,
                     sample_object, sample_timestamp,
                     palette_name, reference_layout):
    '''Create alignment'''
    if int(upload_timestamp) > int(sample_timestamp):
        seq_object = upload_object
        if seq_object is None:
            return ''
        seq_lines = parse_seq_object(seq_object)
    else:
        if sample_object is None:
            return ''
        seq_lines = sample_object

    names, seqs, conservation  = parse_sequences(seq_lines)
 
    x, y, n_seqs, sequence_length = get_dimensions(seqs)

    try:
        ordered_names, ordered_seqs = get_msa_order(reference_name,
                                                    names,
                                                    seqs)
    except KeyError:
        ordered_names, ordered_seqs = names, seqs
    
    palette = COLOR_DIC[palette_name]
    
    # I am sort of misusing the checkbox for the alignment layout. Really this
    # should be returning True/False rather than [True] and []
    
    if reference_layout: 
        reference_layout = True
    else: 
        refrence_layout = False

    text_values, text_colors, block_values, block_colors = \
    alignment_layout(ordered_seqs, layout, palette, reference_layout, BASE_DIC)
    
    trace = go.Heatmap(z=block_values,
                       colorscale = block_colors,
                       showscale=False,
                      )

    steps = [{'args': ['xaxis', {'range': [-0.5 + e, 30.5 + e]}],
              'method': 'relayout',
              'label': ''} for e in range(sequence_length-30)]

    webgl_text = {'type': 'scattergl',
                        'mode': 'text',
                        'x': x,
                        'y': y,
                        'text': text_values,
                        'yaxis': 'y2',
                        'textfont': {
                            'size': 18,
                            'color': text_colors
                        }}

    
    
    bar_trace = {'type': 'bar', 
                        'x': list(range(sequence_length)),
                        'y': conservation,
                        'marker': {'color': '#1e90ff'}}


    fig = tools.make_subplots(rows=2, cols=1,
                              shared_xaxes=True,
                              vertical_spacing=0.001)
    fig.append_trace(trace, 2, 1)
    fig.append_trace(bar_trace, 1, 1)

    fig = fig.to_plotly_json()

    fig['data'].append(webgl_text)

    sliders = [dict(
        minorticklen = 0,
        tickwidth = 0,
        active = 0,
        steps = steps
    )]
    fig['layout'] = dict(
              sliders=sliders,
        yaxis2=dict(autorange='reversed',
                   ticks='',
                   ticksuffix='  ',
                   ticktext=ordered_names,
                   tickvals=list(np.arange(0, len(block_values))),
                   showticklabels=True),
        yaxis=dict(ticks='',
                   ticksuffix='  ',
                   showticklabels=False,
                   domain=[0.7, 1]),
        margin=go.layout.Margin(
            l=200,
            r=50,
            b=0,
            t=50,
            pad=0),
        height=((n_seqs*50) + 100),
        xaxis = {'range': [-0.5, 30.5]},
        showlegend=False
    )
    
    height = (fig['layout']['height'] - fig['layout']['margin']['t'] -
    fig['layout']['margin']['b'])
    y1_height = 65 #px
    fig['layout']['yaxis']['domain'] = [1 - y1_height/height, 1]
    fig['layout']['yaxis2']['domain'] = [0, 1.01 - (y1_height/height)]
    return fig


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', debug=True)
