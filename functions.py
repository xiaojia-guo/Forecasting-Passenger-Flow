#!/usr/bin/env python
__author__ = "Andrey Karasev, Xiaojia Guo"
__license__ = "MIT"
__version__ = "1.1"
__email__ = "andreyjet@gmail.com, x.guo.11@ucl.ac.uk"


import pandas as pd
import numpy as np
import pickle
from sklearn.tree import DecisionTreeRegressor
from scipy.stats import gumbel_r
from bokeh.plotting import figure, show, vplot, hplot, ColumnDataSource, output_file, gridplot
from bokeh.models import HoverTool
from pandas.tseries.offsets import *

pd.options.mode.chained_assignment = None


def dist_generator(leaf, trial_number, coefs):
    """
        Args:
            leaf: terminal node ID 
            trial_number: number of simulations
            coefs: parameters of the Gumbel distribution for each leaf

        Returns:
            an array with shape (1, trial_number) where each value is a random draw from the distribution of the leaf

     """

    g_model = gumbel_r(loc=coefs['loc'][leaf], scale=coefs['scale'][leaf])
    sim_results = g_model.rvs(size=trial_number)

    # Find and replace the non-positive samples.
    for i in range(trial_number):
        if sim_results[i] <= 0:
            sim_results[i] = g_model.rvs(size=1)
    return np.ceil(sim_results)


def stat_calc(quantiles, sim_df, freq, columns):

    """
        Args:
            quantiles: target quantiles
            sim_df: dataFrame with results of simulation (datetime in rows, number of passengers in columns
            freq: forecast resolution

        Returns:
            DataFrame with statistics from simulation, and some varibles needed for making a chart
     """

    indexes = []
    for k, ind in enumerate(sim_df.index):
        if k % freq == 0:
            index = ind
        indexes.append(index)
    df_copy = sim_df.copy()  
    df_copy.index = indexes
    df_group = df_copy.groupby(df_copy.index).sum()
    df_group.index.names = ['inter_beginning']
    # statistics
    quant = quantiles.split(',')
    colnames = ['median','p'+str(quant[0]),'p'+str(quant[1]),'p'+str(quant[2]),'p'+str(quant[3])]
    quant = [float(i) for i in quant]
    colnames = ['median','p'+str(quant[0]),'p'+str(quant[1]),'p'+str(quant[2]),'p'+str(quant[3])]
    df_group[colnames[0]] = df_group[columns].median(axis=1)
    df_group[colnames[1]] = df_group[columns].quantile(q=quant[0]/100, axis=1)
    df_group[colnames[2]] = df_group[columns].quantile(q=quant[1]/100, axis=1)
    df_group[colnames[3]] = df_group[columns].quantile(q=quant[2]/100, axis=1)
    df_group[colnames[4]] = df_group[columns].quantile(q=quant[3]/100, axis=1)
    df_group.drop(columns, axis=1, inplace=True)
    # vars for chart
    df_group['hour'] = df_group.index.hour
    df_group['minute'] = df_group.index.minute
    df_group['inter_end'] = df_group.index + Minute(freq)
    df_group['end_hour'] = df_group['inter_end'].dt.hour
    df_group['end_minute'] = df_group['inter_end'].dt.minute

    return df_group


def data_processing():

    """
    This function grabs and cleans input.xlsx, aircrft_type.csv, coef.csv and treeModel.pickle in the "input" folder.
    It returns the cleaned passenger records and the 
    """
    
    df = pd.read_excel('input/input.xlsx')
    regDict = pd.read_excel('input/regions.xlsx', index_col='IATA').to_dict()
    boyDict = pd.read_csv('input/aircrft_type.csv', index_col='ACRFT_TYPE').to_dict()
    coefs =  pd.read_csv('input/coef.csv', index_col='leaf').to_dict() 
    ukIATA = pd.read_csv('input/UK_IATA.csv')
    with open('input/treeModel.pickle', 'rb') as r:
        model = pickle.load(r)

    # Drop domestic passengers.
    df = df[~df.IB_iata_loc1.isin(ukIATA['UK_IATA'].tolist())]

    # if IB_CHOX_TM is available, script will us this value, if not it will use  IB_EST_CHOX. If last is also
    # unavailable, it will is ATO value
    df['on_chock_best_approx'] = np.nan
    df['on_chock_best_approx'].fillna(df.IB_CHOX_TM, inplace=True)
    df['on_chock_best_approx'].fillna(df.IB_EST_CHOX, inplace=True)
    df['on_chock_best_approx'].fillna(df.ATO, inplace=True)

    # we cannot use cases where we don't have any approximation for on_chock time
    df = df[~df['on_chock_best_approx'].isnull()]

    # first set of features for the tree model
    df['ib_PlanVsOn_chock'] = (df['on_chock_best_approx'] - df['ARR_STO']).astype('timedelta64[m]')
    df['InBoundHour'] = df.on_chock_best_approx.dt.hour
    df['ibFlightLoad'] = df.IB_PAX_TOTAL / df.IB_MAX_PAX
    df['ibFlightLoad'].fillna(1.0, inplace=True)  # in 2% of cases this field is empty
    df['spread'] = (df.OB_STO - df['on_chock_best_approx']).astype('timedelta64[m]')
    df['Ib_Region'] = df.IB_iata_loc1.apply(lambda x: regDict['Retail Markets'][x])
    df['ib_aircraft_body'] = df.ib_arcrft_type.apply(lambda x: boyDict['aircraft_body_type'][x])

    # Transform wide table to long table
    longDf = pd.DataFrame()
    for i, classPax in enumerate(['CLASS1_PAX', 'CLASS2_PAX', 'CLASS3_PAX', 'CLASS4_PAX']):
        temp = df[df[classPax] > 0]
        if len(temp) == 0:
            continue
        globals()['df_' + classPax] = pd.DataFrame()
        for n in range(1, temp[classPax].max() + 1):
            tempClass = temp[temp[classPax] == n]
            globals()['df_' + classPax] = globals()['df_' + classPax].append([tempClass] * n, ignore_index=True)
        globals()['df_' + classPax]['PaxClass'] = i + 1
        longDf = longDf.append(globals()['df_' + classPax])

    # second set of features for the tree model
    longDf['IDAHO_pax_class_3'] = longDf.PaxClass.apply(lambda x: 1 if x == 3 else 0)
    longDf['Ib_Region_N America'] = longDf.Ib_Region.apply(lambda x: 1 if x == 'N America' else 0)
    longDf['ib_aircraft_body_W'] = longDf.ib_aircraft_body.apply(lambda x: 1 if x == 'W' else 0)
    longDf['Ob_Region_Domestic'] = longDf.OB_DICE.apply(lambda x: 1 if x == 'D' else 0)

    inputDf = longDf[['TRANSFER_REC_ID', 'ibFlightLoad', 'InBoundHour', 'spread', 'IDAHO_pax_class_3', 'IB_TERM5_FLAG',
                      'ib_aircraft_body_W', 'IB_STAND_TYPE_P_FLAG', 'Ib_Region_N America', 'Ob_Region_Domestic']]

    passengers = longDf[['TRANSFER_REC_ID','on_chock_best_approx', 'ib_displayed_flyt_no', 'ob_displayed_flyt_no', 'spread']]
    passengers['leaf'] = model.apply(inputDf.drop('TRANSFER_REC_ID'))

    return coefs, passengers


def sim_run(BigStartTime, trialNumber, passengers, columns, coefs, forecast_window, quantiles):

    """
    Function sim_run performs simulation and saves individual results
    to two csv files

        Args:
            BigStartTime: starting time of the forecast window
            trialNumber: number of simulations
            passengers: pre-processed data frame
            columns: list of columns
            coefs: dictionary with parameters of the distribution for each leaf
            forecast_window: size of forecasting window in minutes


        Returns:
           DataFrame with statistics from simulation, and some varibles needed for making a chart.
           Two CSV files - one for each passenger's connection time and another for the expected number
           of late passengers of each outbound flight - will be saved in the "output" folder.
     """
    quant = quantiles.split(',')
    IbUpperLimit = BigStartTime + Minute(forecast_window)
    IbLowerLimit = BigStartTime - Minute(150)
    chartDf = pd.DataFrame(index=pd.DatetimeIndex(start=BigStartTime, end=IbUpperLimit, freq='min'))

    workingDf = passengers[((passengers.on_chock_best_approx >= IbLowerLimit)
                            & (passengers.on_chock_best_approx <= IbUpperLimit))]
    workingDf.reset_index(drop=True, inplace=True)

    deltaMins = list(map(lambda x: dist_generator(x, trialNumber, coefs), workingDf.leaf.tolist()))
    simOut = pd.DataFrame(deltaMins, dtype=int)
    temp = pd.concat([workingDf[['TRANSFER_REC_ID','on_chock_best_approx', 'spread', 'ob_displayed_flyt_no']], simOut], axis=1)
    if workingDf.shape[0] == 0:
        chartDf = pd.DataFrame(np.zeros((chartDf.shape[0], trialNumber)),
                               index=pd.DatetimeIndex(start=BigStartTime, end=IbUpperLimit, freq='min'))
        colnames = ['median','p'+str(quant[0]),'p'+str(quant[1]),'p'+str(quant[2]),'p'+str(quant[3])]

        temp = pd.DataFrame(columns=colnames)
        temp.drop('spread',axis=1).to_csv('output/individual_%s.csv' % str(BigStartTime).replace(" ", "_").replace(":", "-"))

    else:
        late_pass = temp[['spread', 'ob_displayed_flyt_no']].groupby('ob_displayed_flyt_no')[['spread']].count()
        for shift in [0, 5, 10, 15]:
            temp_min = temp.copy()
            for col in columns:# I limited number columns otherwise it take to long to run
                temp_min[col] = temp_min[[col,'spread']].apply(lambda x: 1 if x[col] > (x['spread'] - 35 + shift) else 0, axis=1)
            group_temp_min = temp_min.drop(['TRANSFER_REC_ID','on_chock_best_approx', 'spread'], axis=1).groupby('ob_displayed_flyt_no').sum()
            late_pass['+'+str(shift)+' min'] = group_temp_min[columns].median(axis=1)
        late_pass.drop('spread',axis=1).to_csv('output/late_passengers_%s.csv' % str(BigStartTime).replace(" ", "_").replace(":", "-"))

        for col in columns:
            temp[col] = pd.to_timedelta(temp[col], 'm') + temp.on_chock_best_approx
            groupedTrial = pd.DataFrame(temp[col].value_counts())
            groupedTrialR = groupedTrial.reindex(
            index=pd.DatetimeIndex(start=BigStartTime, end=IbUpperLimit, freq='min'), fill_value=0)
            chartDf = chartDf.merge(groupedTrialR, left_index=True, right_index=True)
        quant = quantiles.split(',')
        colnames = ['median','p'+str(quant[0]),'p'+str(quant[1]),'p'+str(quant[2]),'p'+str(quant[3])]
        quant = [float(i) for i in quant]
        
        temp[columns] = temp[columns].astype('int64')
        temp[colnames[0]] = temp[columns].median(axis=1)
        temp[colnames[1]] = temp[columns].quantile(q=quant[0]/100, axis=1)
        temp[colnames[2]] = temp[columns].quantile(q=quant[1]/100, axis=1)
        temp[colnames[3]] = temp[columns].quantile(q=quant[2]/100, axis=1)
        temp[colnames[4]] = temp[columns].quantile(q=quant[3]/100, axis=1)
        prob = np.zeros(simOut.shape[0])
        for i in range(0, simOut.shape[0]):
            thisRow = simOut.loc[i,:]
            prob[i] = round(len(thisRow[thisRow > (workingDf['spread'][i] - 35)])/trialNumber, 2)
        temp['P(missing connecting flight)'] = prob
        temp.drop(columns, axis=1, inplace=True)
        temp['ib_displayed_flyt_no'] = workingDf.ib_displayed_flyt_no
        temp['ob_displayed_flyt_no'] = workingDf.ob_displayed_flyt_no
    for col in colnames[0:5]:
        temp[col] = pd.to_datetime(temp[col])
        temp.drop('spread',1).to_csv('output/individual_%s.csv' % str(BigStartTime).replace(" ", "_").replace(":", "-"))
    return chartDf


def results(sim_results, columns, StartTime, quantiles, resolutions=[1, 5, 15, 60]):
    
    """
        Function results provides the aggregated predictions -the transfer passenger flow at the RtF desk.
        It saves the results in seperate CSV files
        
        Args:
        sim_results: data frame with results of simulation (date/time in rows, number of passengers in columns
        resolutions: list with time resolutions
        columns: list of columns
        StartTime: starting time of a forecast window
        
        Returns:
        CSV files and some varibles that are useful for generating the passenger flow figures.
        """
    
    for r, res in enumerate(resolutions):
        globals()['chartDf%smin' % res]=stat_calc(quantiles, sim_results, res, columns)
        quant = quantiles.split(',')
        quant = [float(i) for i in quant]
        colnames = ['median','p'+str(quant[0]),'p'+str(quant[1]),'p'+str(quant[2]),'p'+str(quant[3])]
        colnames.insert(0,'inter_end')
        globals()['chartDf%smin' % res][colnames].to_csv('output/outPut_%smin_%s.csv' % (res,str(StartTime).replace(" ", "_").replace(":", "-")))
    
    return [globals()['chartDf%smin' % res] for res in resolutions]

def chart(dfs, StartTime, resolutions, quantiles):

    """
    Function generates html dashboard

        Args:
            dfs: data frames with simulation statistics
            StartTime: starting time of the forecast window
            resolutions: list with time resolutions

     """

    TOOLS = "pan,reset,hover,save"
    quant = quantiles.split(',')
    quant = [float(i) for i in quant]
    colnames = ['median','p'+str(quant[0]),'p'+str(quant[1]),'p'+str(quant[2]),'p'+str(quant[3])]
   
    for i, data in enumerate(dfs):
        source = ColumnDataSource(data=dict(
            median=data['median'],
            p05=data[colnames[1]],
            p95=data[colnames[4]],
            left=data.index,
            right=data.inter_end,
            hour=data['hour'],
            minute=data['minute'],
            hour_end=data['end_hour'],
            minute_end=data['end_minute']
            ))
        globals()['p'+str(i)] = figure(tools=TOOLS, plot_width=600, plot_height=400,
                                       title='%s min. intervals' % resolutions[i], x_axis_type="datetime")
        globals()['p'+str(i)].title_text_font_size = '20pt'
        globals()['p'+str(i)].yaxis.axis_label = "Number of passengers"
        globals()['p'+str(i)].xaxis.axis_label = "Time"
        quant = quantiles.split(',')
        cols = ['median','p'+str(quant[3]),'p'+str(quant[0])]
        globals()['p'+str(i)].quad(top='p95', bottom=0, left='left', right='right', source=source,color="#dadaeb",
                                   alpha=0.80)
        globals()['p'+str(i)].quad(top='median', bottom=0, left='left', right='right', source=source,color="#6a51a3",
                                   fill_alpha=0.80)
        globals()['p'+str(i)].quad(top='p05', bottom=0, left='left', right='right', source=source, color="#3f007d",
                                   fill_alpha=0.80)
        globals()['p'+str(i)].xgrid.grid_line_color = None
        globals()['p'+str(i)].ygrid.grid_line_alpha = 0.5
        globals()['p'+str(i)].ygrid.grid_line_dash = [6, 4]

        hover = globals()['p'+str(i)].select_one(HoverTool)
        hover.point_policy = "follow_mouse"
        hover.tooltips = [ ("Time interval", " @hour:@minute - @hour_end:@minute_end "),
                          (quant[0] + "% chance to have less than", " @p95{‘0’} pass."),
                          ("50% chance to have at least", " @median{‘0’} pass."),
                           (quant[3] + "% chance to have less than", " @p05{‘0’} pass.")]

    plot_list = [globals()['p'+str(plot)] for plot in range(len(dfs))]
    plot_list_trans = []
    for plot_number, plot in enumerate(plot_list):
        if plot_number % 2 == 0:
            try:
                plot_list_trans.append([plot_list[plot_number], plot_list[plot_number+1]])
            except:
                plot_list_trans.append([plot_list[plot_number]])
        else:
            continue

    grid=gridplot(plot_list_trans)

    output_file("output/plot%s.html" % str(StartTime).replace(" ", "_").replace(":", "-"))
    show(grid)
