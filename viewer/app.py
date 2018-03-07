import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import skbio
from skbio.sequence import DNA
import numpy as np
import base64
import io
import datetime
import plotly.figure_factory as ff
import plotly.graph_objs as go
from collections import OrderedDict
import json

#modify this path for different alignment file. This will be replaced in the
#the future with the uploasd component
test_data_fp = '../data/msa1.fna'
msa = skbio.alignment.TabularMSA.read(test_data_fp, constructor=DNA)
base_dic = {'A': 1, 'C': .25, 'G': .5, 'T': .75, '-': 0}

def get_tick_vals(start, stop):
    tick_vals = np.arange(start, stop, 10)
    tick_vals = np.array([9, 19, 29]) - (tick_vals%10)
    return tick_vals


names = [seq.metadata['id'] for seq in msa]
colorscale_blocks = [[0.00, '#F4F0E4'],
              [0.25, '#1b9e77'],
              [0.50, '#d95f02'],
              [0.75, '#7570b3'],
              [1.00, '#e7298a']]

colorscale_letters = [[0.00, '#FFF'], 
                    [0.25, '#FFF'], 
                    [0.50, '#FFF'], 
                    [0.75, '#FFF'],
                    [1.00, '#FFF']]
letter_colors = {'A': '#e7298a', 'C': '#1b9e77', 'G': '#d95f02', 'T': '#7570b3', '-': '#444'}

app = dash.Dash()

app.layout = html.Div(children=[
    html.H1(children='Alignment Viewer'),
    html.Div([ 
    dcc.RadioItems(
                id='layout-type',
                options=[{'label': i, 'value': i} for i in ['Blocks', 'Letters']],
                value='Blocks',
                labelStyle={'display': 'inline-block'}
            ),

    html.Div([
    html.Label('Reference Sequence'),
    dcc.Dropdown(
        id='parent-seq',
        options=[{'label': label, 'value': label} for label in names],
        value=names[0])], style={'width': '22%', 'display': 'inline-block'}),

    dcc.Graph(id='alignment',
        config={
            'displayModeBar': False}
        ),
    dcc.Slider(
        id='alignment-slider',
        updatemode='drag',
        min=0,
        max=len(msa[0]) - 30,
        value=0,
        step=1
        )
    ], style={'width': '80%', 'display': 'inline-block'}),
    html.Div(id='ordered-values', style={'display': 'none'})
    ]
)


@app.callback(
        dash.dependencies.Output('ordered-values', 'children'), 
        [dash.dependencies.Input('parent-seq', 'value')])
def seq_align_for_plot(name):
    seq_dic = OrderedDict(zip(names, msa))
    seq_dic.move_to_end(name, last=False)
    base_text = [list(str(seq_dic[e])) for e in seq_dic]
    base_values = np.zeros((len(base_text), len(base_text[0])))
    for i in range(len(base_text[0])):
        for j in range(len(base_text)):
            if base_text[j][i] != base_text[0][i]:
                base_values[j][i] = base_dic[base_text[j][i]]
    updated_names = list(seq_dic.keys())
    jdump = json.dumps([base_text, base_values.tolist(), updated_names])
    return(jdump)


@app.callback(
    dash.dependencies.Output('alignment', 'figure'),
    [dash.dependencies.Input('alignment-slider', 'value'),
     dash.dependencies.Input('layout-type', 'value'),
     dash.dependencies.Input('ordered-values', 'children')])
def update_figure(start, layout, bases):
    stop = start + 30
    tick_values = get_tick_vals(start, stop)
    if layout == 'Blocks':
        colorscale = colorscale_blocks
        font_properties = {
                'family': 'Courier New, monospace',
                'size': 14,
                'color': '#3f566d'}
    else:
        colorscale = colorscale_letters
        font_properties = {
                'family': 'Courier New Bold, monospace',
                'size':16,
                'color': '#000'
                }
    base_text, base_values, ordered_names = json.loads(bases)
    base_subset = np.array(base_values)[:, start:stop]
    text_subset = np.array(base_text)[:, start:stop]
    
    fig = ff.create_annotated_heatmap(base_subset,
                                      y=ordered_names,
                                      annotation_text=text_subset,
                                      colorscale=colorscale)
    
    fig['layout'].update(
            xaxis=dict(side='top',       
                       ticktext=np.arange(start, stop, 10) - (np.arange(start, stop, 10) % 10) + 10,
                       tickvals=tick_values,
                       showticklabels=True,
                       tickfont=dict(family='Bookman',
                                     size=18,
                                     color='#22293B',
                                    ),
                       ),
            
            yaxis=dict(autorange='reversed',
               ticks='',
               ticksuffix='  ',
               ticktext=ordered_names,
               tickvals=list(np.arange(0, len(base_text))),
               showticklabels=True),
            
        width=1200,
        height=(25*len(base_subset)),
        margin=go.Margin(
            l=200,
            r=50,
            b=0,
            t=50,
            pad=0),

        annotations=dict(font=font_properties)
        )
    if layout == 'Letters':
        for letter in fig.layout.annotations:
            letter['font']['color'] = letter_colors[letter['text']]
            if letter['text'] != '-' and letter['y'] != ordered_names[0]:
                if letter['text'] == text_subset[0][letter['x']]:
                    letter['text'] = '.'
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
