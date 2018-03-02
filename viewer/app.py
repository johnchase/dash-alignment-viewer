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


msa = skbio.alignment.TabularMSA.read('/Users/jc33/Dropbox/plotly/msa10.fna', constructor=DNA)
base_dic = {'A': 1, 'C': .25, 'G': .5, 'T': .75, '-': 0}

def seq_align_for_plot(msa):
    base_text = [list(str(e)) for e in msa]
    base_values = np.zeros((len(base_text), len(base_text[0])))
    for i in range(len(base_text[0])):
        for j in range(len(base_text)):
            if base_text[j][i] != base_text[0][i]:
                base_values[j][i] = base_dic[base_text[j][i]]
    return(base_text, base_values)

def get_tick_vals(start, stop):
    tick_vals = np.arange(start, stop, 10)
    tick_vals = np.array([9, 19, 29]) - (tick_vals%10)
    return tick_vals

base_text, base_values = seq_align_for_plot(msa)
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
    dcc.Dropdown(
        id='parent-seq',
        options=[{'label': label, 'value': label} for label in names],
        value=names[0])], style={'width': '25%', 'display': 'inline-block'}),

    dcc.Graph(id='alignment',
        config={
            'displayModeBar': False}
        ),
    dcc.Slider(
        id='alignment-slider',
        updatemode='drag',
        min=0,
        max=len(base_values[0]) - 30,
        value=0,
        step=1
        )
    ], style={'width': '80%', 'display': 'inline-block'})
    ]
)

@app.callback(
    dash.dependencies.Output('alignment', 'figure'),
    [dash.dependencies.Input('alignment-slider', 'value'),
     dash.dependencies.Input('layout-type', 'value')])
def update_figure(start, layout):
    stop = start + 30
    tick_values = get_tick_vals(start, stop)
    print(layout) 
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

    base_subset = np.array(base_values)[:, start:stop]
    text_subset = np.array(base_text)[:, start:stop]
    
    fig = ff.create_annotated_heatmap(base_subset,
                                      y=names,
                                      annotation_text=text_subset,
                                      colorscale=colorscale)
    
    fig['layout'].update(
            xaxis=dict(ticks='',
                       side='top',       
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
               ticktext=names,
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
            if letter['text'] != '-' and letter['y'] != '1AM1JR7QWMSFA_8958':
                if letter['text'] == text_subset[0][letter['x']]:
                    letter['text'] = '.'
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
