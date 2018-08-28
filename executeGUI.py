#!/usr/bin/env python
__author__ = "Andrey Karasev"
__license__ = "MIT"
__version__ = "1.0"
__email__ = "andreyjet@gmail.com"


import pandas as pd
import functions as fc
import time
from apscheduler.schedulers.blocking import BlockingScheduler  # need to install with pip
from pandas.tseries.offsets import *
from tkinter import *
import tkinter as tk
from tkinter import ttk


sched = BlockingScheduler()
counter = 0


def get_values():
    forecast_window = forecast_windowVar.get() - 1 # to make window more intuitive. Use can put 120 min.
    trialNumber = trialNumberVar.get()
    resolutions = resolutionsVar.get().split(',')
    resolutions = [int(number) for number in resolutions]
    start = startVar.get()
    end = endVar.get()
    frequency = frequencyVar.get()
    quantiles = quantilesVar.get()

    @sched.scheduled_job('interval', minutes=frequency, start_date=start, end_date=end)
    def main():
        t0 = time.time()
        columns = list(range(trialNumber))
        timePoint = pd.DatetimeIndex(freq='%smin' % frequency, start=start, end=end)
        global counter
        StartTime = timePoint[counter]# + Hour(4)  # Day(2)
        print('Working on output for %s' % str(StartTime))
        coefs, passengers = fc.data_processing()
        chartDf = fc.sim_run(StartTime, trialNumber, passengers, columns, coefs, forecast_window, quantiles)
        dfList = fc.results(chartDf, columns, StartTime, quantiles, resolutions)
        fc.chart(dfList, StartTime, resolutions, quantiles)
        t1 = time.time()
        counter = counter + 1
        print('Done with outputs for %s in %s seconds' % (str(StartTime), int(t1 - t0)))

    root.destroy()

root = Tk()
root.wm_title("Heathrow_app")

forecast_windowVar = IntVar()
forecast_windowVar.set(120)

quantilesVar = StringVar()
quantilesVar.set('5,10,90,95')

trialNumberVar = IntVar()
trialNumberVar.set(500)

resolutionsVar = StringVar()
resolutionsVar.set('5,15,60')

frequencyVar = IntVar()
frequencyVar.set(15)

startVar = StringVar()
startVar.set(str(pd.Timestamp.now()).split('.')[0])

endVar = StringVar()
endVar.set(str(pd.Timestamp.now()+Day()).split('.')[0])


for i in range(7):
    globals()['spacer'+str(i)] =ttk.Label(root, text="            ").grid(column=i, row=0)

for i in range(12):
    globals()['spacerR'+str(i)] =ttk.Label(root, text=" ").grid(column=0, row=i)


s = ttk.Style()
s.configure('.', font=('Helvetica', 16))
title = ttk.Label(root, text="Parameters:            ", font=("Helvetica", 20)).grid(column=1, row=1, columnspan=1)
title1 = ttk.Label(root, text="Forecasting window (min)                                              ",
                  font=("Helvetica", 14)).grid(column=1, row=3, columnspan=4)
title2 = ttk.Label(root, text="Number of simulations                                                 ",
                  font=("Helvetica", 14)).grid(column=1, row=4, columnspan=4)
title3 = ttk.Label(root, text="Resolutions (numbers separated with commas)            ",
                   font=("Helvetica", 14)).grid(column=1, row=5, columnspan=4)
title4 = ttk.Label(root, text="Quantiles (numbers separated with commas)                  ",
                   font=("Helvetica", 14)).grid(column=1, row=6, columnspan=4)
title5 = ttk.Label(root, text="Update frequency (min)                                                 ",
                  font=("Helvetica", 14)).grid(column=1, row=7, columnspan=4)
title6 = ttk.Label(root, text="Starting time (YYYY-MM-DD HH:MM:SS)                        ",
                  font=("Helvetica", 14)).grid(column=1, row=8, columnspan=4)
title7 = ttk.Label(root, text="Ending time (YYYY-MM-DD HH:MM:SS)                          ",
                  font=("Helvetica", 14)).grid( column=1, row=9, columnspan=4)


F_W = ttk.Entry(root, textvariable=forecast_windowVar, font=("Helvetica", 14)).grid(column=5, row=3, columnspan=1)
T_N = ttk.Entry(root, textvariable=trialNumberVar, font=("Helvetica", 14)).grid(column=5, row=4, columnspan=1)
RS = ttk.Entry(root, textvariable=resolutionsVar, font=("Helvetica", 14)).grid(column=5, row=5, columnspan=1)
QUANTILES = ttk.Entry(root, textvariable=quantilesVar, font=("Helvetica", 14)).grid(column=5, row=6, columnspan=1)
FQ = ttk.Entry(root, textvariable=frequencyVar, font=("Helvetica", 14)).grid(column=5, row=7, columnspan=1)
START = ttk.Entry(root, textvariable=startVar, font=("Helvetica", 14)).grid(column=5, row=8, columnspan=1)
END = ttk.Entry(root, textvariable=endVar, font=("Helvetica", 14)).grid(column=5, row=9, columnspan=1)

RUN = tk.Button(root, text='START', relief='raised', command=get_values, font=("Helvetica", 20)).grid(
    column=5, row=10, columnspan=1, rowspan=1)

root.mainloop()
print('Script is started at %s' % pd.Timestamp.now())
sched.start()
