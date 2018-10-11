"""Functions for application
"""
import collections
import base64
import skbio

import pandas as pd

from io import StringIO


def parse_sequences(seq_lines): 
    msa = skbio.alignment.TabularMSA.read(seq_lines,
                                          constructor=skbio.sequence.DNA)
    seqs, names = zip(*[(str(seq), seq.metadata['id']) for seq in msa])
    conservation = msa.conservation()
    return names, seqs, conservation


def alignment_layout(seqs, layout_type, letter_colors, base_dic):
    '''Get layout for alignment'''
    letter_colors['-'] = '#444'
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

                        [0.1, letter_colors['C']],
                        [0.26, letter_colors['C']],

                        [0.26, letter_colors['G']],
                        [0.51, letter_colors['G']],

                        [0.51, letter_colors['T']],
                        [0.76, letter_colors['T']],

                        [0.76, letter_colors['A']],
                        [1.00, letter_colors['A']]]

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
    try:
        _, content_string = seq_object.split(',')
        decoded = base64.b64decode(content_string)
        seq_lines = StringIO(decoded.decode('utf-8'))
        return seq_lines
    except AttributeError:
        pass
