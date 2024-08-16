import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

def plot_many_genes(df, genes, s=25, figwidth=4, figheight=8, numcols=5, title=None, dpi=None, show_zero_wells=False, title_fontsize=20):
    
    n = len(genes)
    
    num_rows = ( (n-1)//numcols ) + 1
    num_cols = np.min([n,numcols])
    
    df = np.array(df.loc[:,["x_position", "y_position"] + genes])
    
    if not dpi==None:
        plt.subplots(num_rows, num_cols, figsize=(figwidth*num_cols, figheight*num_rows), dpi=dpi)
    else:
        plt.subplots(num_rows, num_cols, figsize=(figwidth*num_cols, figheight*num_rows))
    
    for i in range(num_rows * num_cols):
        
        plt.subplot(num_rows, num_cols, i+1)
        
        # If gene[i], make the plot
        if i < n:
            if not show_zero_wells:
                sample = df[df[:,i+2]>0]
            
            sns.scatterplot(x=sample[:,0], y=sample[:,1], hue=sample[:,i+2], s=s)
            plt.axis("equal")
            plt.title(genes[i])
        
        # We want to remove ticks and frame regardless
        plt.xticks([])
        plt.yticks([])
        plt.box(False)
    
    if not title==None:
        plt.suptitle(title, fontsize=title_fontsize)
        
    plt.show()