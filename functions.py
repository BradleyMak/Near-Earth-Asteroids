#import relevant modules
from astroquery.jplhorizons import Horizons
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import functions
import subprocess
import os
import time

def read_from_excel(object, obs_code = 995):
    """Reads astrometric data from Excel file with the name '{Object} Data' and saves it as a txt file with the correct format for find_orb.

    Args:
        object (string): Object name as given in the Excel filename.
        obs_code (int, optional): Observatory code. Defaults to 995 (Durham). Optional as it is not needed if Excel spreadsheet contains observatory information.
    
    Returns:
        pandas dataframe: Dataframe containing Excel data.
    """
    #reading data from excel file
    df = pd.read_excel(f'{object} Data.xlsx')
    #print(df)

    #creating (or overwriting) new file with a relevant name.
    f = open(f"{object}_data.txt", "w")
    f.writelines(['123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|\n', '<-ObjDesig->*nnYYYY MM DD.DDDDD HH MM SS.SSSsdd mm ss.ss<blanks >MM.MMBz<ref>COD\n'])

    #lots of string slicing to take data from the excel spreadsheet and write it to the txt file line by line in the correct format.
    object_list = []
    for j in range(0, 12):
        try:
            object_list.append(object[j])
        except:
            object_list.append(' ')
    object_designation = ''.join(map(str, object_list))

    for i in range(0, len(df)):
        Date_and_Time = str(df['Date + Time'][i])
        year = Date_and_Time[0:4]
        if Date_and_Time == 'nan':
            continue
        month = Date_and_Time[5:7]; day = Date_and_Time[8:10]+'.'+str(np.round((float(Date_and_Time[11:13])*3600 + float(Date_and_Time[14:16])*60 + float(Date_and_Time[17:len(Date_and_Time)]))/(24*3600), 5))[2:7]
        Date_and_Time = f'{year} {month} {day}'

        RA = df['RA'][i]
        hours = RA[0:2]; minutes = RA[3:5]; seconds = str(np.round(float(RA[6:len(RA)]), 2))
        while seconds[2] != '.':
            seconds = '0' + seconds
        while len(seconds) != 5:
            seconds += '0'
        RA_err = str(np.round(float(df['RA err (as)'][i]), 3))
        while len(RA_err) != 5:
            RA_err += '0'

        Dec = df['Dec'][i]
        if Dec[0] != '+' and Dec[0] != '-':
            degrees = '+'+Dec[0:2]; arcmins = Dec[3:5]; arcsecs = Dec[6:len(Dec)]
        else:
            degrees = Dec[0:3]; arcmins = Dec[4:6]; arcsecs = Dec[7:len(Dec)]
        Dec_err = str(np.round(float(df['Dec err (as)'][i]), 3))
        while len(Dec_err) != 5:
            Dec_err += '0'

        #if Excel file contains column for observatory code, read the observatory code...
        try:
            obs_code = df['Observatory'][i]
        #... if not, default to using the observatory code given as an argument to the function.
        except:
            obs_code = obs_code    

        f.write(f'{object_designation}   {Date_and_Time} {hours} {minutes} {seconds} {degrees} {arcmins} {arcsecs}                     {obs_code} {RA_err} {Dec_err}\n')

    #closes and saves file.
    f.close()

    return df

def dec_to_dms(decimal):
    """Converts Dec in degrees to sexigesimal format.

    Args:
        decimal (string): Takes the Dec in deg as a string as a parameter.

    Returns:
        string: Returns Dec in format '+/-hh mm ss' as a string.
    """
    #checks whether Dec is positive or negative.
    if decimal < 0:
        sign = '-'
    else:
        sign = '+'
    decimal = abs(decimal)
    degrees = int(decimal // 1)
    minutes = int(((decimal % 1)*60) // 1)
    seconds = np.round((((decimal % 1)*60) % 1)*60, 2)
    #if statements to return Dec in correct format.
    if degrees < 10:
        degrees = f'0{degrees}'
    if minutes < 10:
        minutes = f'0{minutes}'
    if seconds < 10:
        seconds = f'0{seconds}'
    return f'{sign}{degrees} {minutes} {seconds}'

def dec_to_deg(Dec):
    """Converts declination from sexigesimal to degrees.

    Args:
        Dec (string): Declination in format dd:mm:ss.sssss...

    Returns:
        float: Declination in degrees.
    """
    degrees = Dec[0:2]; minutes = Dec[3:5]; seconds = Dec[6:len(Dec)]
    return float(degrees)+float(minutes)/60+float(seconds)/3600

def RA_to_deg(RA):
    """Converts right ascension from sexigesimal to degrees.

    Args:
        Dec (string): Right ascension in format hh:mm:ss.sssss...

    Returns:
        float: Right ascension in degrees.
    """
    hours = RA[0:2]; minutes = RA[3:5]; seconds = RA[6:len(RA)]
    return float(hours)*15+15*float(minutes)/60+15*float(seconds)/3600

def RA_to_hms(decimal):
    """Converts RA in degrees to sexigesimal format.

    Args:
        decimal (string): Takes the RA in deg as a string as a parameter.

    Returns:
        string: Returns RA in format 'hh mm ss' as a string.
    """
    hours = int(decimal // 15)
    minutes = int(((decimal % 15)*60) // 15)
    seconds = np.round(((((decimal % 15)*60) % 15)*60)/15, 2)
    #if statements to return RA in correct format.
    if hours < 10:
        hours = f'0{hours}'
    if minutes < 10:
        minutes = f'0{minutes}'
    if seconds < 10:
        seconds = f'0{seconds}'
    return f'{hours} {minutes} {seconds}'

def month_to_number(month):
    """Converts month as a word into its corresponding integer (01-12).

    Args:
        month (string): Month (abberviation)

    Returns:
        string: Month (integer enclosed within a string)
    """
    #converts month in writing to its corresponding integer.
    return {
            'jan': '01',
            'feb': '02',
            'mar': '03',
            'apr': '04',
            'may': '05',
            'jun': '06',
            'jul': '07',
            'aug': '08',
            'sep': '09', 
            'oct': '10',
            'nov': '11',
            'dec': '12'
    }[month.lower()]

def run_find_orb(observations_filename):
    """Runs find_orb software.

    Args:
        observations_filename (string): Filename of the file containing the astrometry to be passed into find_orb.
    """
    #sets current working directory to find_c64 directory. This will need to be changed for other users to the relevant directory containing the fo64 executable along with all of the other find_orb files.
    os.chdir('C:\\Users\\bradl\\find_c64')
    #runs find_orb executable with given observations file.
    subprocess.run([".\\fo64.exe", f'C:\\Users\\bradl\\.vscode\\advanced_lab\\{observations_filename}'], shell=True)

def read_fo_elements(elements_filename):
    """Function to read orbital elements from a txt file in the format that find_orb returns.

    Args:
        elements_filename (string): Filename of the file containing the orbital elements in the format as returned by find_orb.

    Returns:
        strings: Orbital elements returned by find_orb.
    """
    f = open(f'{elements_filename}')
    lines = f.readlines()
    f.close()
    fo_peri_arg = str.split(lines[4])[5]
    fo_sma = str.split(lines[5])[1]
    fo_asc_node = str.split(lines[5])[5]
    fo_eccentricity = str.split(lines[6])[1]
    fo_inc = str.split(lines[6])[5]
    fo_period = str.split(lines[7])[1]
    fo_peri_dist = str.split(lines[8])[1]
    fo_apogee_dist = str.split(lines[8])[5]
    return fo_peri_arg, fo_sma, fo_asc_node, fo_eccentricity, fo_inc, fo_period, fo_peri_dist, fo_apogee_dist

def get_jpl_ephemeris(object_designation, obs_code, start_date, end_date, step):
    """Gets JPL ephemeris for a given object and range of dates.

    Args:
        object_designation (string): Object name/designation as given in JPL Horizons.
        obs_code (string/integer): Observatory code.
        start_date (string): Ephemeris start date given in format 'YYYY-MM-DD'.
        end_date (string): Ephemeris end date given in format 'YYYY-MM-DD'.
        step (string): Time between ephemeris entries given in format e.g. '3d' for 3 days.

    Returns:
        HorizonsEphemeris: JPL Horizons ephemeris for requested date range.
    """
    eph = Horizons(id=f'{object_designation}', location = f'{obs_code}', epochs = {'start':f'{start_date}', 'stop':f'{end_date}', 'step': f'{step}'}).ephemerides(quantities = 1)
    return eph

def get_jpl_elements(object_designation, location = '500@10'):
    """Reads orbital elements (relative to the Sun, by default) for a given object from JPL Horizons.

    Args:
        object_designation (string): Object name/designation as given in JPL Horizons.
        location (string, optional): Location the orbital parameters are given relative to. Defaults to '500@10' (the Sun).

    Returns:
        string: Orbital elements returned by JPL Horizons.
    """
    elements = Horizons(id=f'{object_designation}', location = location).elements()
    #print(elements.columns)
    jpl_eccentricity = elements['e'][0]
    jpl_sma = elements['a'][0]
    jpl_period = str(float(elements['P'][0])/365.25)
    jpl_peri_dist = elements['q'][0]
    jpl_apogee_dist = elements['Q'][0]
    jpl_peri_arg = elements['w'][0]
    jpl_asc_node = elements['Omega'][0]
    jpl_inc = elements['incl'][0]
    return jpl_peri_arg, jpl_sma, jpl_asc_node, jpl_eccentricity, jpl_inc, jpl_period, jpl_peri_dist, jpl_apogee_dist

def read_jpl_ephemeris(eph, object_designation, obs_code, rows_to_read = 0):
    """Reads a given JPL ephemeris and converts it into the correct format for find_orb and saves it to the file 'jpl_eph.txt'.

    Args:
        eph (HorizonsEphemeris): JPL Horizons ephemeris.
        object_designation (string): Object name/designation as given in JPL Horizons.
        obs_code (string/integer): Observatory code.
        rows_to_read (int, optional): Number of rows of the ephemeris the user would like to be read and written to the resulting txt file. Defaults to 0, which tells the function to read all rows.
    """
    #if user has not specified a certain number of ephemeris rows to read into the txt file, then by default read all rows.
    if rows_to_read == 0:
        rows_to_read = len(eph)
    #sets current working directory to be the the advanced_lab directory. This will need to be changed for other users to the relevant project directory.
    os.chdir('C:\\Users\\bradl\\.vscode\\advanced_lab')
    f = open("jpl_eph.txt", "w")
    f.writelines(['123456789|123456789|123456789|123456789|123456789|123456789|123456789|123456789|\n', '<-ObjDesig->*nnYYYY MM DD.DDDDD HH MM SS.SSSsdd mm ss.ss<blanks >MM.MMBz<ref>COD\n'])

    #iterates over all rows that the user has requested to be read.
    for i in range(0, rows_to_read):
        #converts date and time into correct format.
        jpl_date = eph['datetime_str'][i]
        year = jpl_date[0:4]; month = jpl_date[5:8]; day = jpl_date[9:11] + '.' + str(np.round((int(jpl_date[12:14])*3600+int(jpl_date[15:17])*3600)/(24*3600), 5))[2:]
        if len(day) < 8:
            day = day.ljust(8, '0')
        print(jpl_date)
        date = f'{year} {month_to_number(month)} {day}'

        #Adds zeroes to the left of the RA and Dec to give them the correct length (for the find_orb file format).
        jpl_RA = RA_to_hms(eph['RA'][i]).ljust(11, '0')
        jpl_DEC = dec_to_dms(eph['DEC'][i]).ljust(12, '0')

        #writes line to jpl ephemeris txt file in the correct format (hence the excessive whitespace).
        f.write(f'{object_designation}   {date} {jpl_RA} {jpl_DEC}                     {obs_code}\n')

    #closes and saves file.
    f.close()