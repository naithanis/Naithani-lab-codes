#This script was written by Daemon Dikeman and revised +tested by Sushma Naithani, Oregon State University

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys


sys.argv = ['diffX_cluster.py']
#XL_namebase = r"C://Users//Naithans//Desktop//filename.csv"
metrix = 'euclidean'
if len(sys.argv) > 1:
    if sys.argv[1] == '-m':
        metrix = sys.argv[2]
        if len(sys.argv) > 3:
            XL_Name = sys.argv[3]
    else:
        XL_Name = sys.argv[1]

XL_handle = open(XL_Name,  "r")
#df_Diff_XL = ["","","",""]
df_Diff_XL = pd.read_csv(XL_handle,
                            index_col=0,
                            na_values="",
                            )
XL_handle.close()
# df_Diff_XL.drop(labels=["GeneIDs", "Design Element"], axis=1, inplace=True)


df_Diff_XL.fillna(value=0, inplace=True)

name = "Filename"
curfig = [6, 6] #figsize
metrica=metrix
xticks=["B1", "B2", 'B3', 'B4'] 

Diff_XL_clus = sns.clustermap(df_Diff_XL,
                                metric=metrica,
                                cmap='RdBu_r', #this color schema, where RdBu is red to blue
                                figsize=curfig,
                                col_cluster=False,
                                xticklabels=True,
                                yticklabels=True,
                                cbar_kws={"label": 'Log2 Fold-Change',
                                        'orientation': 'horizontal'}, #for bar it could be vertical
                                center=0.0,
                                vmax=5.0
                                
                                # vmax=6.0 : 
)


top, bottom = Diff_XL_clus.ax_heatmap.get_ylim()
Diff_XL_clus.ax_heatmap.set_ylim(top + 0.5, bottom - 0.5)
Diff_XL_clus.ax_col_dendrogram.set_visible(False)
BaseDendroBox = Diff_XL_clus.ax_col_dendrogram.get_position()
BaseDendroBox.y0 = 1.01 * BaseDendroBox.y0
BaseDendroBox.y1 = (BaseDendroBox.y1 + 8 * BaseDendroBox.y0) / 9
Diff_XL_clus.cax.set_position (BaseDendroBox)
Diff_XL_clus.ax_heatmap.set_ylabel("Gene IDs")
Diff_XL_clus.ax_heatmap.tick_params(axis='y', labelsize=8) 
Diff_XL_clus.ax_heatmap.tick_params(axis='x', labelsize=10) 
Diff_XL_clus.cax.xaxis.set_ticks_position("top")
Diff_XL_clus.cax.xaxis.set_label_position("top")
plt.savefig("/Users/downloads/FigureX.png",
            bbox_inches='tight'
            )
