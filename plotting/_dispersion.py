import _plotter

class Dispersion(_plotter.plotter):
    def __init__(self,
            dir1,
            dir2,
            title1='Before',
            title2='After',
            filename='dispersion',
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
        '''
        Plot the dispersion plot
        '''
        _plotter.plotter.__init__(self,dir1,dir2,macro='dispersion.gp')
        self._set_ranges(xmin,xmax,ymin,ymax,ymin2,ymax2)
        self.title1=title1
        self.title2=title2

        self.energy_in_gev=energy

        self._init_filename(out_folder,filename)
        if plot:
            self.plot()
