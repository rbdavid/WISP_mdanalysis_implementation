
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from sklearn.cluster import KMeans

# ----------------------------------------
# PLOTTING FUNCTIONS:
# ----------------------------------------

def plot_square_matrix(square_matrix,figure_name,axes_label='',cbar_label='',plotting_cmap='bwr',v_range=None,minor_ticks=10,major_ticks=100):
        """
        """
        nNodes = len(square_matrix)
        node_range = range(nNodes+1)
        fig, ax = plt.subplots()
        ax.tick_params(which='major',length=6,width=2)
        ax.tick_params(which='minor',length=3,width=1)
        ax.xaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.xaxis.set_major_locator(MultipleLocator(major_ticks))
        ax.yaxis.set_minor_locator(MultipleLocator(minor_ticks))
        ax.yaxis.set_major_locator(MultipleLocator(major_ticks))
        
        if v_range != None:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap,vmin=v_range[0],vmax=v_range[1])
        else:
                temp = plt.pcolormesh(node_range,node_range,square_matrix,cmap=plotting_cmap)
        cb1 = plt.colorbar()
        cb1.set_label(r'%s'%(cbar_label))

        xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:]]
        ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:]]
        temp.axes.set_xticks(temp.axes.get_xticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_xticks(temp.axes.get_xticks()[:]+0.5)
        temp.axes.set_yticks(temp.axes.get_yticks(minor=True)[:]+0.5,minor=True)
        temp.axes.set_yticks(temp.axes.get_yticks()[:]+0.5)
        temp.axes.set_xticklabels(xlabels)
        temp.axes.set_yticklabels(ylabels)

        plt.xlim((-0.5,nNodes+0.5))
        plt.ylim((-0.5,nNodes+0.5))
        plt.xlabel(axes_label,size=14)
        plt.ylabel(axes_label,size=14)
        ax.set_aspect('equal')
        plt.tight_layout()
        plt.savefig(figure_name,dpi=600,transparent=True)
        plt.close()

def paths_analysis_plotting(paths,source_indices,sink_indices,nNodes,length_freq_file_name,node_freq_file_name,bins=100): 
        """
        """
        lengths = [path[0] for path in paths]
        plt.hist(lengths,bins = bins)
        plt.xlabel('Path Lengths',size=14)
        plt.ylabel('Frequency',size=14)
        plt.xlim((0,np.max(lengths)))
        plt.tight_layout()
        plt.savefig(length_freq_file_name,dpi=600,transparent=True)
        plt.close()

        node_histogram = np.zeros((nNodes),dtype=np.int)
        for path in paths:
                for node in path[1:]:
                        node_histogram[node] += 1

        fig, ax = plt.subplots()
        plt.bar(range(nNodes),node_histogram,width=1.0)
        plt.xlabel('Node Index',size=14)
        plt.ylabel('Frequency',size=14)

        for source in source_indices:
                ax.annotate('Source\nNode %d'%(source),verticalalignment='bottom',horizontalalignment='left',xy=(source+0.5,node_histogram[source]),xytext=(source+10.5,node_histogram[source]-10.0),arrowprops=dict(width=5.0,headwidth=0.0,facecolor='black', shrink=0.05))       # headlength=0.0,

        for sink in sink_indices:
                ax.annotate('Sink\nNode %d'%(sink),verticalalignment='bottom',horizontalalignment='left',xy=(sink+0.5,node_histogram[sink]),xytext=(sink+10.5,node_histogram[sink]-10.0),arrowprops=dict(width=5.0,headwidth=0.0,facecolor='black', shrink=0.05))       # headlength=0.0,

        plt.xlim((-0.5,nNodes))
        plt.ylim((0,len(paths)*1.10))
        plt.tight_layout()
        plt.savefig(node_freq_file_name,dpi=600,transparent=True)
        plt.close()


