"""App style parameters
"""
import colorlover as cl

cl_scales = cl.scales['4']['qual']
COLOR_DIC = {}
for palette in cl_scales.keys():
    COLOR_DIC[palette] = {key: value for key, value 
                          in zip(list('ATGC'), cl_scales[palette])}


BASE_DIC = {'A': 1, 'C': .25, 'G': .5, 'T': .75, '-': 0}

UPLOAD_BUTTON = {'width': '100',
                 'height': '30',
                 'background': 'transparent',
                 'border': '2px solid #0099CC',
                 'margin': '2px',
                 'border-radius': '6px',
                 'outline': 'none',
                 'text-align': 'center',
                 'vertical-align': 'middle',
                 'display': 'inline-block'
               }
