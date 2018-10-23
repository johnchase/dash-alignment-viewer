"""App style parameters
"""
import colorlover as cl

cl_scales = cl.scales['4']['qual']
COLOR_DIC = {}


COLOR_DIC['DRuMS Nucleic Acid'] = {'A': 'rgb(80,80,255)', 'C': 'rgb(224,0,0)',
        'T': 'rgb(230,230,0)', 'G': 'rgb(0,192,0)'}

COLOR_DIC['Jalview'] = {'A': 'rgb(100,247,63)', 'C': 'rgb(255,180,48)',
        'G': 'rgb(235,65,60)', 'T': 'rgb(60,136,238)'}

COLOR_DIC['Purines/Pyrimidine'] = {'A': 'rgb(255,131,250)', 
        'C': 'rgb(64,224,208)', 'G': 'rgb(255,131,250)', 
        'T': 'rgb(64,224,208)'}

for palette in cl_scales.keys():
    COLOR_DIC[palette] = {key: value for key, value 
                          in zip(list('ATGC'), cl_scales[palette])}


BASE_DIC = {'A': 1, 'C': .25, 'G': .5, 'T': .75, '-': 0}

MENU_ELEMENTS = {'width': '150px', 
                 'marginRight': '20px', 
                 'display': 'inline-block'}

TITLES = {'fontWeight': 'bold', 'fontSize': 20, 'color': '#4C4C4C'} 

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


