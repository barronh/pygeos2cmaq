__all__ = 'simpledate minus1hour'.split()

from datetime import datetime, timedelta

def simpledate(date, p):
    return date.strftime(p)

def minus1hour(date, p):
    return (date - timedelta(hours = 1)).strftime(p)
