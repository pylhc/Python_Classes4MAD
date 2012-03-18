import os



macros_path=os.path.join(os.path.dirname(__file__),'macros')

class plotter:
    def __init__(self,
            dir1,
            dir2,
            macro,
            beam=1):
        self.dir1=dir1
        self.dir2=dir2
        self.macro=macro
        self.fname_end='.eps'
        self.terminal='postscript enhanced color solid 22'
        self.beam=beam
        self.xtics_plot2=self._get_xticks()


    def __getitem__(self,value):
        '''
        So that the plotter works as
        a dictionary..
        '''
        return getattr(self,value)

    def plot(self,filename=None,tmp_file='gplot.tmp.gp'):
        '''
        Function plots to filename
        usually called by init functions
        of child classes.
        '''
        if filename:
            self.filename=filename
        script=file(tmp_file,'w')
        mac=file(self.get_macro_file(),'r').read()
        script.write(mac % self)
        script.close()
        os.system('gnuplot '+tmp_file)

    def get_macro_file(self):
        '''
        Returns the full path to
        the macro file
        '''
        return os.path.join(macros_path,self.macro) 

    def _set_ranges(self,xmin,xmax,ymin,ymax,ymin2=None,ymax2=None):
        '''
        Sets the ranges. This function is meant to
        be used internally by the init functions
        '''
        self.xmin=0.0
        self.xmax=26600
        self.ymin=ymin
        self.ymax=ymax

        if ymin2==None:
            self.ymin2=ymin
        else:
            self.ymin2=ymin2

        if ymax2==None:
            self.ymax2=ymax
        else:
            self.ymax2=ymax2

        if self.ymax=='*' or self.ymin2=='*':
            return 0
        if self.ymax==-self.ymin2:
            self.ymax-=0.01

    def _init_filename(self,out_folder,filename):
        '''
         Used by the __init__ funtions
        '''
        if out_folder:
            if not os.path.isdir(out_folder):
                os.makedirs(out_folder)
            self.filename=os.path.join(out_folder,filename+self.fname_end)
        else:
            self.filename=filename+self.fname_end
    def _get_xticks(self):
        fpath=os.path.join(macros_path,'xtics_b'+str(self.beam)+'.gp')
        return file(fpath,'r').read()
