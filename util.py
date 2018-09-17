"""Functions for application
"""
import collections
import base64


import pandas as pd


def parse_sequences(seq_lines):
    """Parse sequences that are present as lines i.e. from reading a file"""
    names = []
    seqs = []
    for seq_line in seq_lines:
        if seq_line.startswith('>'):
            names.append(seq_line.split(' ')[0].replace('>', ''))
        else:
            seqs.append(seq_line.strip())
    return names, seqs


def alignment_layout(seqs, layout_type, letter_colors, base_dic):
    '''Get layout for alignment'''
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
                        [0.1, '#F4F0E4'],

                        [0.1, '#1b9e77'],
                        [0.26, '#1b9e77'],

                        [0.26, '#d95f02'],
                        [0.51, '#d95f02'],

                        [0.51, '#7570b3'],
                        [0.76, '#7570b3'],

                        [0.76, '#e7298a'],
                        [1.00, '#e7298a']]

    for seq in seqs[1:]:
        seq_series = pd.Series(list(seq))

        if layout_type == 'Letter':
            text_colors.extend(seq_series.replace(letter_colors).tolist())
            seq_series.where(seq_series != pd.Series(list(seqs[0])), '.',
                             inplace=True)
            text_values.extend(seq_series.tolist())

        elif layout_type == 'Block':
            seq_series.where(seq_series != pd.Series(list(seqs[0])), 0,
                             inplace=True)
            seq_series.replace(base_dic, inplace=True)
            block_values.append(seq_series.tolist())
    return text_values, text_colors, block_values, block_colors


def get_msa_order(reference_name, names, seqs):
    """Order sequences"""
    seq_dic = collections.OrderedDict(zip(names, seqs))
    seq_dic.move_to_end(reference_name)
    return zip(*list(seq_dic.items())[::-1])

def get_dimensions(seqs):
    '''Function for getting figure get_dimensions'''
    n_seqs, sequence_length = len(seqs), len(seqs[0])
    y = [[i]*sequence_length for i in range(n_seqs)]
    y = [item for sublist in y for item in sublist]
    x = list(range(sequence_length))*n_seqs
    return x, y, n_seqs, sequence_length

def parse_seq_object(seq_object):
    '''Parse the object that is read in by Dash. This would ultimately use
    skbio hoewever that causes Heroku to fail currently'''
    _, content_string = seq_object.split(',')
    decoded = base64.b64decode(content_string)
    seq_lines = decoded.decode('utf-8').strip().split('\n')
    return parse_sequences(seq_lines)
