## Exceptions

__all__ = ['PyDSTool_BoundsError', 'PyDSTool_KeyError',
           'PyDSTool_UncertainValueError', 'PyDSTool_TypeError',
           'PyDSTool_ExistError', 'PyDSTool_AttributeError',
           'PyDSTool_ValueError', 'PyDSTool_UndefinedError',
           'PyDSTool_InitError']


class PyDSTool_Errors(Exception):
    pass


class PyDSTool_UncertainValueError(PyDSTool_Errors):
    def __init__(self, value, varval=None):
        if varval is None:
            valstr = ''
        else:
            valstr = ' at variable = '+str(varval)
        self.varval = varval
        self.value = value+valstr;
    def __str__(self):
        return repr(self.value)

class PyDSTool_BoundsError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)

class PyDSTool_KeyError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)

class PyDSTool_ValueError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)

class PyDSTool_TypeError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)
    
class PyDSTool_ExistError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)

class PyDSTool_UndefinedError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)

class PyDSTool_AttributeError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)

class PyDSTool_InitError(PyDSTool_Errors):
    def __init__(self, value):
        self.value = value;
    def __str__(self):
        return repr(self.value)
