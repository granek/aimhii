import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter, MaxNLocator

import os
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

READ_HEIGHT = 1
READ_HEIGHT_INCHES = 0.07
READ_PITCH = READ_HEIGHT * 1.3
MARGIN=1
def plot_clusters(metacluster_list,plotfilename):

    plotbase,plotext = os.path.splitext(plotfilename)
    for i, metacluster in enumerate(metacluster_list):
        logger.info("drawing metacluster: [0]".format(metacluster))
        rect_list,stack_count = metacluster.draw()
        fig = plt.figure(figsize=(8, READ_HEIGHT_INCHES*stack_count+MARGIN))
        ax = plt.axes(frameon=False)
        ax.get_xaxis().tick_bottom() # Turn off ticks at top of plot
        ax.axes.get_yaxis().set_visible(False) # Hide y axis


        # ax = fig.add_subplot(111)
        for rect in rect_list:
            ax.add_patch(rect)
        ax.set_xlim(metacluster.iv.start,metacluster.iv.end)
        ax.set_ylim(1,stack_count*READ_PITCH)
        ax.set_xlabel(metacluster.iv.chrom)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f'%x))
        xmin, xmax = ax.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
        plt.tight_layout()
        plt.savefig("{0}_{1}{2}".format(plotbase,i,plotext))
