import _plotter

class NormDispersion(_plotter.plotter):
    '''
    Plot the normalized dispersion.

    Only dir1 and dir2 are required input.

    :param dir1: Folder where result 1 is located
    :param dir2: Folder where result 2 is located
    :param title1: Data set title for result 1
    :param title2: Data set title for result 2
    :param filename: filename of plot (no file ending)
    :param out_folder: folder where plot should be placed
    :param xmin: x-axis minimum value
    :param xmax: x-axis maximum value
    :param ymin: y-axis minimum value
    :param ymax: y-axis maximum value
    :param ymin2: y-axis minimum value, 2. plot
    :param ymax2: y-axis maximum value, 2. plot
    '''
    def __init__(self,
            dir1,
            dir2,
            title1='Before',
            title2='After',
            filename='normalized_dispersion',
            out_folder='',
            xmin=0.0,
            xmax=26600.0,
            energy=450,
            ymin='*',
            ymax='*',
            ymin2=None,
            ymax2=None,
            plot=True
            ):

        _plotter.plotter.__init__(self,dir1,dir2,macro='normdispersion.gp')

        self._set_ranges(xmin,xmax,ymin,ymax,ymin2,ymax2)

        self.title1=title1
        self.title2=title2

        self.energy_in_gev=energy

        self._init_filename(out_folder,filename)

        if plot:
            self.plot()
