'''
postRunClean removes unnecessary files after an NLL run, and generates the
final NLL result file. It also generates a map of the events.
'''

import sys
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)

from dataclasses import dataclass
import pygmt as pg
import pandas as pd
from complem_figures.basemap_lines import ABOVE_SEA_COLOR, PYGMT_SHORELINE, draw_relief_pygmt

@dataclass
class MapPostRunParams:
    fileBulletin: str
    figSave: str

# FUNCTION
def removeHighErr(df):
    df = df[df.erh <= 3.0]
    df = df[df.erv <= 3.0]
    df = df[df.gap <= 300]
    df = df[df.rms <= 0.5]
    return df

def genFigure(parameters):
    #---- Read OBS file
    with open(parameters.fileBulletin,'r',encoding='utf-8') as f:
        lines = f.readlines()
    print(f"Catalog succesfully read @ {parameters.fileBulletin}")

    events = [event.split() for event in lines]
    events_df = pd.DataFrame(events).drop(columns=[0,1,2,3,4,5,11]).rename(columns={6:'Latitude',7:'Longitude',8:'Depth',9:'Magnitude',10:'rms',12:'erh',13:'erv',14:'gap'}).astype(float)
    events_df = removeHighErr(events_df)

    #---- Set Pyrenees borders
    region = [-2.25,3.5,42,44]

    fig = pg.Figure()
    with pg.config(MAP_FRAME_TYPE="fancy+"):
        fig.basemap(region=region, projection="M6i", frame='af')
    draw_relief_pygmt(fig, region)
    fig.coast(shorelines=PYGMT_SHORELINE, resolution='i', area_thresh='0/0/1', borders="1/0.75p,black")

    #---- Plot events
    with pg.config(COLOR_BACKGROUND='black', COLOR_FOREGROUND=ABOVE_SEA_COLOR):  # swapped by reverse
        pg.makecpt(cmap="viridis", series=[0,15,1], reverse=True)
    fig.plot(
        x=events_df.Longitude,
        y=events_df.Latitude,
        style="c0.02c",
        # size=1 * events_df.Magnitude, # the Magnitude is 0 for all
        fill=events_df.Depth,
        cmap=True,
        transparency=30,
    )

    fig.colorbar(frame=['a5f5+lDepth [km] (deeper than 15 in black, above sea level in magenta)'])
    fig.savefig(parameters.figSave, dpi=300)

    print(f"Figure succesfully saved @ {parameters.figSave}")


# MAIN
if __name__ == "__main__":
    params_W_40 = MapPostRunParams(
        fileBulletin = 'RESULT/GLOBAL_W_40.txt',
        figSave = 'RESULT/MAPS/GLOBAL_W_40.pdf',
    )

    genFigure(params_W_40)

    params_W_80 = MapPostRunParams(
        fileBulletin = 'RESULT/GLOBAL_W_80.txt',
        figSave = 'RESULT/MAPS/GLOBAL_W_80.pdf',
    )

    genFigure(params_W_80)

    params_W_140 = MapPostRunParams(
        fileBulletin = 'RESULT/GLOBAL_W_140.txt',
        figSave = 'RESULT/MAPS/GLOBAL_W_140.pdf',
    )

    genFigure(params_W_140)

    params_W_200 = MapPostRunParams(
        fileBulletin = 'RESULT/GLOBAL_W_200.txt',
        figSave = 'RESULT/MAPS/GLOBAL_W_200.pdf',
    )

    genFigure(params_W_200)