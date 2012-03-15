import _plotter


class betabeat(_plotter.plotter):
    def __init__(self,
            dir1,
            dir2,
            title1='Before',
            title2='After',
            filename='betabeat',
            xmin=0.0,
            xmax=26600.0,
            ymin='*',
            ymax='*',
            y2min=None,
            y2max=None,
            plot=True
            ):
        '''
        Plot the famous betabeat plot
        '''
        super(betabeat,self).__init__(dir1,dir2,macro='betabeat.gp')
        self._set_ranges(xmin,xmax,ymin,ymax,y2min,y2max)
        if plot:
            self.plot(filename+self.fname_end)
        else:
            self.filename=filename+self.fname_end
