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
            title='LHCB1 450 GeV',
            title1='Before',
            title2='After',
            beam=1,
            filename='normalized_dispersion',
            out_folder='',
            xmin=0.0,
            xmax=26600.0,
            ymin='*',
            ymax='*',
            ymin2=None,
            ymax2=None,
            free_1=False,
            free_2=False,
            plot=True
            ):

        if dir2:
            macro='base.gp'
        else:
            macro='base1.gp'

        _plotter.plotter.__init__(self,dir1,dir2,macro=macro)

        self._set_ranges(xmin,xmax,ymin,ymax,ymin2,ymax2)

        self.title=title
        self.title1=title1
        self.title2=title2

        self._set_yminl()

        self._set_datasets(free_1,free_2,'ND','D')

        self._init_filename(out_folder,filename)

        self.ylabel1="{/Symbol D}D_{x}/{/Symbol b}_x^{1/2} [m^{1/2}]"
        self.ylabel2="{/Symbol D}D_{y} [m]"

        self.xfunc="2"
        self.yfunc="($4-$8)"
        self.yfunc2="($4-$7)"
        self.errfunc="5"

        if plot:
            self.plot()
