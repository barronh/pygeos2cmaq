from yaml import load as yload
from datetime import datetime, timedelta

def loader(config = 'config.yaml'):
    """
    Process configuration dictionary to modify types
    """
    
    # Load YAML
    config_dict = yload(file(config))
    
    # Convert dates using YYYY-MM-DD HH:MM:SS
    exec('start_date = datetime.strptime(start_date, "%Y-%m-%d %H:%M:%S")', None, config_dict)
    exec('end_date = datetime.strptime(end_date, "%Y-%m-%d %H:%M:%S")', None, config_dict)
    
    # Convert time increment assuming that the unit is a keyword
    config_dict.setdefault('time_incr', '1 hours')
    exec('time_incr = timedelta(%s = %s)' % tuple(config_dict['time_incr'].split(' ')[::-1]), None, config_dict)
    config_dict.setdefault('repair_aero', True)
    config_dict.setdefault('zero_negs', True)
    config_dict.setdefault('no_clobber', False)
    return config_dict

if __name__ == '__main__':
    x = loader()