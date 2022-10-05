# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 18:41:22 2022

@author: carlo
"""
import pygrib
import numpy as np
import pandas as pd
from functools import reduce

#grib = 'dati_copernicus.grib' # Set the file name of your input GRIB file
grib = 'Copernicus_data_h11.grib'
grbs = pygrib.open(grib)      # open the file


# obtain list of all entries
grb_list = grbs.select()
n_list = len(grb_list)

# list of all attributes
grb_list[0].keys()

# find latitude and longitude of a single point of the list
grb_list[3].latlons()


# variable names
var_names = list(set([grb.parameterName for grb in grb_list]))




# support lists for each info I have
df_dewpoint = []
df_temp = []
df_skin = []
df_solar = []


# list of df for each variable
list_df = [[] for _ in var_names]

for grb in grb_list:
    curr_name = grb.parameterName
    line = {
        'year': grb.year,                              # metti anno
        'month': grb.month,                            # metti mese
        'day': grb.day,                                # metti giorno
        'hour': grb.validityTime,                      # metti ora
        'value': grb.values                            # metti valore
        }
    list_df[var_names.index(curr_name)].append(line)


list_df = [pd.DataFrame(df) for df in list_df]    # convert each list to dataframe


# update names of variables
for i in range(len(var_names)):
    list_df[i] = list_df[i].rename(columns={'value':var_names[i]})


# merge by 'year', 'month', 'day', 'hour' the dataframes in a unique one
df_all = reduce(lambda left, right:     # Merge DataFrames in list
                     pd.merge(left , right,
                              on = ['year', 'month', 'day', 'hour'],
                              how = 'inner'),
                     list_df)


# save in csv format
df_all.to_csv('dati_copernicus_09_05_22.csv', sep=',', header=True, index=False)


 