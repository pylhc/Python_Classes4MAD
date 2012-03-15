import os

macros_path=os.path.join(os.path.dirname(__file__),'macros')

class plotter:
    def __init__(self,
            dir1,
            dir2,
            macro):
        self.dir1=dir1
        self.dir2=dir2
        self.macro=macro
        self.fname_end='.eps'
        self.terminal='postscript enhanced color solid 22'

    def __getitem__(self,value):
        '''
        So that the plotter works as
        a dictionary..
        '''
        return getattr(self,value)

    def plot(filename,tmp_file='gplot.tmp'):
        '''
        Function plots to filename
        usually called by init functions
        of child classes.
        '''
        self.filename=filename
        script=file(tmp_file,'w')
        mac=file(self.get_macro_file(),'r')
        script.write(mac % self)
        script.close()
        os.system('gnuplot '+tmp_name)

    def _set_ranges(xmin,xmax,ymin,ymax,y2min=None,y2max=None):
        '''
        Sets the ranges. This function is meant to
        be used internally by the init functions
        '''
        self.xmin=0.0
        self.xmax=26600
        if y2min==None:
            self.y2min=ymin
        else:
            self.y2min=y2min
        if y2max==None:
            self.y2max=ymax
        else:
            self.y2max=y2max

    def get_macro_file():
        '''
        Returns the full path to
        the macro file
        '''
        return os.path.join(macros_path,self.macro) 



