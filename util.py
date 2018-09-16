import pandas as pd

def alignment_layout(seqs, layout_type, letter_colors, base_dic):
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
