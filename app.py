import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
# import skbio
# from skbio.sequence import DNA
import pandas as pd
from collections import OrderedDict
import numpy as np


test_data_fp = 'data/msa10.fna'
lines = tuple(open(test_data_fp, 'r'))
names = []
seqs = []
for line in lines:
    if line.startswith('>'):
        names.append(line.split(' ')[0].replace('>', ''))
    else:
        seqs.append(line.strip())
# msa = skbio.alignment.TabularMSA.read(test_data_fp, constructor=DNA)
# names = [seq.metadata['id'] for seq in msa]


n_seqs, sequence_length = len(seqs), len(seqs[0])
y = [[i]*sequence_length for i in range(n_seqs)]
y = [item for sublist in y for item in sublist]
x = list(range(sequence_length))*n_seqs

letter_colors = {'A': '#e7298a', 'C': '#1b9e77', 'G': '#d95f02', 'T': '#7570b3', '-': '#444'}
base_dic = {'A': 1, 'C': .25, 'G': .5, 'T': .75, '-': 0}

def alignment_layout(seqs, layout_type):
    text_values = list(seqs[0])
    block_values = [[0]*len(text_values)]


    if layout_type == 'Letter':
        text_colors = pd.Series(text_values).replace(letter_colors).tolist()
        block_colors = [[0.00, '#FFF'],
                        [1.00, '#FFF']]
        block_values *= len(seqs)

    elif layout_type == 'Block':
        text_colors = '#000'
        text_values = list(''.join(seqs))
        block_colors = [[0.00, '#F4F0E4'],
                        [0.1,  '#F4F0E4'],

                        [0.1, '#1b9e77'],
                        [0.26, '#1b9e77'],

                        [0.26, '#d95f02'],
                        [0.51, '#d95f02'],

                        [0.51, '#7570b3'],
                        [0.76, '#7570b3'],

                        [0.76, '#e7298a'],
                        [1.00, '#e7298a']]

    for seq in seqs[1:]:
        s = pd.Series(list(seq))

        if layout_type == 'Letter':
            text_colors.extend(s.replace(letter_colors).tolist())
            s.where(s != pd.Series(list(seqs[0])), '.', inplace=True)
            text_values.extend(s.tolist())

        elif layout_type == 'Block':
            s.where(s != pd.Series(list(seqs[0])), 0, inplace=True)
            s.replace(base_dic, inplace=True)
            block_values.append(s.tolist())
    return text_values, text_colors, block_values, block_colors


def get_msa_order(reference_name, names, seqs):
    seq_dic = OrderedDict(zip(names, seqs))
    seq_dic.move_to_end(reference_name)
    return zip(*list(seq_dic.items())[::-1])

app = dash.Dash()

app.layout = html.Div(children=[
    html.H1(children='Alignment Viewer'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),
    html.Div([
    dcc.RadioItems(
    id='layout-type',
    options=[{'label': i, 'value': i} for i in ['Block', 'Letter']],
    value='Block',
    labelStyle={'display': 'inline-block'}),

    html.Div([
    html.Label('Reference Sequence'),
    dcc.Dropdown(
        id='parent-seq',
        options=[{'label': label, 'value': label} for label in names],
        value=names[0])], style={'width': '22%', 'display': 'inline-block'}),

    dcc.Graph(
        id='alignment',
        config={
            'displayModeBar': False}
    )
])])


@app.callback(
    dash.dependencies.Output('alignment', 'figure'),
    [dash.dependencies.Input('layout-type', 'value'),
     dash.dependencies.Input('parent-seq', 'value')])
def create_alignment(layout, reference_name):

    ordered_names, ordered_seqs = get_msa_order(reference_name, names, seqs)
    text_values, text_colors, block_values, block_colors = alignment_layout(ordered_seqs, layout)

    trace = go.Heatmap(z=block_values,
                       colorscale = block_colors,
                       showscale=False,
                      )

    steps = [{'args': ['xaxis', {'range': [-0.5 + e, 30.5 + e]}],
              'method': 'relayout',
              'label': ''} for e in range(sequence_length-30)]

    data=[trace]

    data.append({'type': 'scattergl',
                        'mode': 'text',
                        'x': x,
                        'y': y,
                        'text': text_values,
                        'textfont': {
                            'size': 14,
                            'color': text_colors,
                            'family': 'Helvetica'
                        }})

    sliders = [dict(
        minorticklen = 0,
        tickwidth = 0,
        active = 0,
        steps = steps
    )]

    layout = dict( sliders=sliders,
    yaxis=dict(autorange='reversed',
                   ticks='',
                   ticksuffix='  ',
                   ticktext=ordered_names,
                   tickvals=list(np.arange(0, len(block_values))),
                   showticklabels=True),
        margin=go.layout.Margin(
            l=200,
            r=50,
            b=0,
            t=50,
            pad=0),
        height=(n_seqs*50),
        xaxis = {'range': [-0.5, 30.5]}
    )

    fig = dict(data=data, layout=layout)
    return fig



if __name__ == '__main__':
    app.run_server(debug=False)
