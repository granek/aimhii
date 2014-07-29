import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

READ_HEIGHT = 1
READ_PITCH = READ_HEIGHT * 1.2
def plot_clusters(metacluster_list,plotfilename):

    plotbase,plotext = os.path.splitext(plotfilename)
    for i, metacluster in enumerate(metacluster_list):
        plt.figure()
        ax = plt.gca()
        print "drawing metacluster:", metacluster
        rect_list,stack_count = metacluster.draw()
        for rect in rect_list:
            ax.add_patch(rect)
        ax.set_xlim(metacluster.iv.start,metacluster.iv.end)
        ax.set_ylim(1,stack_count*READ_PITCH)
        plt.savefig("{0}_{1}{2}".format(plotbase,i,plotext))
