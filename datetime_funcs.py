# -*- coding: utf-8 -*-
# @Author: Martin Grunnill
# @Date:   2023-10-03 11:25:43
# @Last Modified by:   Martin Grunnill
# @Last Modified time: 2023-10-03 16:45:12

from dateutil.parser import parse
from datetime import datetime as dt
from datetime import timedelta
import time

def is_date(string, fuzzy=False):
    """
    Return whether the string can be interpreted as a date.

    string: str
        String to check for date.
    fuzzy: bool
        Ignore unknown tokens in string if True
    """
    try: 
        parse(string, fuzzy=fuzzy)
        return True

    except ValueError:
        return False


def date_to_year_fraction(date):
    """Convert datetime.datetime object to year fraction.

    Parameters
    ----------
    date : datetime.datetime
        
    """
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def str_date_to_year_fraction(date_string, delimiter='-'):
    """Convert string date object to year fraction.

    Assumes YYYY to day format

    Parameters
    ----------
    date_string : datetime.datetime
        
    """
    if len(date_string)==4:
        if not date_string.isdigit():
            raise ValueError('If date_string (' + date_string + ') is 4 characters they must be digits')
        return date_string + '.5'
    else:
        if date_string.count(delimiter)==2:
            date_string_split = date_string.split(delimiter)
            if (len(date_string_split[0])!=4 or not date_string_split[0].isdigit() or
                    len(date_string_split[1])!=2 or not date_string_split[1].isdigit()):
                raise ValueError('String_date (' + date_string + ') is not recognised as a YYYY'+delimiter+'MM format.')
            date_string = date_string + delimiter + '15'

        try:
            date_string = parse(date_string, yearfirst=True)
        except:
            raise ValueError('String_Date (' + date_string + ') is not recognised as a date format.')

        return date_to_year_fraction(date_string)

def string_convert_date_format(string, old_format='%m-%d-%Y', new_format='%Y-%m-%d'):
    return dt.strptime(string, old_format).strftime(new_format)


def year_fraction_to_date(year_fraction):
    """ Convert year fraction to datetime object.

        Parameters
        ----------
        year_fractiom : float

        Returns
        -------
        datetime
    """
    year = int(year_fraction)
    remainder = year_fraction - year
    boy = dt(year, 1, 1)
    eoy = dt(year + 1, 1, 1)
    seconds = remainder * (eoy - boy).total_seconds()
    return boy + timedelta(seconds=seconds)