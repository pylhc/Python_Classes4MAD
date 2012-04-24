
## 
# This deals with the 
# append functionality
# of the script writing..

__append__=True
def __get_append__():
    return __append__

def set_append(append):
    '''
     Set to True (default) to 
     make one long plot script
     containing all plots
    ''' 
    if append not in [True,False]:
        raise ValueError('append must be True or False')
    global __append__
    __append__=append
