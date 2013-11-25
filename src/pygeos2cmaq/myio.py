__all__ = ['myio']
import warnings

def formatwarning(message, category, filename, lineno, line = 0):
    strout = "\n***********\n%s:%s: %s:\n\t%s\n***********\n" % (filename, lineno, category.__name__, message)
    return strout

warnings.formatwarning = formatwarning

class myio(object):
    def __init__(self):
        self._warnings = ""
        self._errors = ""
        self._status = ""
    def warn(self, *args, **kwds):
        keep = kwds.pop('keep', True)
        if 'message' in kwds:
            message = kwds['message']
        else:
            message = args[0]
        if keep:
            self._warnings += '\n' + message
        warnings.warn(*args, **kwds)
    
    def status(self, *args, **kwds):
        show = kwds.get('show', True)
        if 'message' in kwds:
            message = kwds['message']
        else:
            message = args[0]
        self._status += '\n' + message
        if show: print message

    def error(self, *args, **kwds):
        if 'message' in kwds:
            message = kwds['message']
        else:
            message = args[0]
        self._errors += '\n' + message
        self.warn(message, keep = False)
    
    def getwarnings(self):
        return self._warnings
    
    def geterrors(self):
        return self._errors
    
    def getstatus(self):
        return self._status