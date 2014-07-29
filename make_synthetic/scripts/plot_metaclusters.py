import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import os

READ_HEIGHT = 1
READ_HEIGHT_INCHES = 0.08
READ_PITCH = READ_HEIGHT * 1.2
MARGIN=1
def plot_clusters(metacluster_list,plotfilename):

    plotbase,plotext = os.path.splitext(plotfilename)
    for i, metacluster in enumerate(metacluster_list):
        print "drawing metacluster:", metacluster
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
        xmin, xmax = ax.get_xaxis().get_view_interval()
        ymin, ymax = ax.get_yaxis().get_view_interval()
        ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
        plt.tight_layout()
        plt.savefig("{0}_{1}{2}".format(plotbase,i,plotext))

        #============================================================
        # fig1 = plt.figure(facecolor='white')
        # ax1 = plt.axes(frameon=False)

        # #ax1.set_frame_on(False) # Alternate way to turn frame off

        # ax1.get_xaxis().tick_bottom() # Turn off ticks at top of plot

        # #ax1.axes.get_xaxis().set_visible(False)
        # ax1.axes.get_yaxis().set_visible(False) # Hide y axis

        # # Add a plot
        # y_offset = 2.0
        # x = numpy.arange(-5.0, 5.0, 0.1)
        # for i in range(3):
        #     ax1.plot(x, numpy.sin((i + 1) * x) + i * y_offset, label=str(i+1))

        # plt.legend(loc='lower right')

        # # Draw the x axis line
        # # Note that this must be done after plotting, to get the correct
        # # view interval
        # xmin, xmax = ax1.get_xaxis().get_view_interval()
        # ymin, ymax = ax1.get_yaxis().get_view_interval()
        # ax1.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))
        # #============================================================



        
