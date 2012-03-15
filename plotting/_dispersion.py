import _plotter

class dispersion(_plotter.plotter):
    def __init__(self,
            dir1,
            dir2,
            title1='Before',
            title2='After',
            filename='dispersion',
            xmin=0.0,
            xmax=26600.0,
            ymin='*',
            ymax='*',
            y2min=None,
            y2max=None
            ):
        '''
        Plot the dispersion plot
        '''
        super(betabeat,self).__init__(dir1,dir2,macro='dispersion.gp')
        self._set_ranges(xmin,xmax,ymin,ymax,y2min,y2max)
        self.filename=filename+self.fname_end
        self.plot(filename+self.fname_end)
