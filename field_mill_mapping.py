'''
All data from https://kscwxarchive.ksc.nasa.gov/
'''
author = "@rtphokie"
from mpl_toolkits.basemap import Basemap
from pprint import pprint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import unittest
import simple_cache
from matplotlib.pyplot import figure
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from basemapcache import BasemapWrapper
import os
from datetime import datetime, date, time
from tqdm import tqdm
from pytz import timezone
from matplotlib import gridspec
import matplotlib.dates as mdates

threshold = 1500
tz = timezone('US/Eastern')

day = 'wed'
if day == 'wed':
    t0 = datetime(2020, 5, 27, 16, 33)
    limit = datetime(2020, 5, 27, 16, 18)
if day == 'sat':
    t0 = datetime(2020, 5, 30, 15, 22)
    limit = datetime(2020, 5, 30, 15, 7)

t0 = tz.localize(t0)
limit = tz.localize(limit)
print (t0)


def calc_new_coord(lat1, lon1, rad, dist):
    """
    Calculate coordinate pair given starting point, radial and distance
    Method from: http://www.geomidpoint.com/destination/calculation.html
    """

    flat = 298.257223563
    a = 2 * 6378137.00
    b = 2 * 6356752.3142

    # Calculate the destination point using Vincenty's formula
    f = 1 / flat
    sb = np.sin(rad)
    cb = np.cos(rad)
    tu1 = (1 - f) * np.tan(lat1)
    cu1 = 1 / np.sqrt((1 + tu1*tu1))
    su1 = tu1 * cu1
    s2 = np.arctan2(tu1, cb)
    sa = cu1 * sb
    csa = 1 - sa * sa
    us = csa * (a * a - b * b) / (b * b)
    A = 1 + us / 16384 * (4096 + us * (-768 + us * (320 - 175 * us)))
    B = us / 1024 * (256 + us * (-128 + us * (74 - 47 * us)))
    s1 = dist / (b * A)
    s1p = 2 * np.pi

    while (abs(s1 - s1p) > 1e-12):
        cs1m = np.cos(2 * s2 + s1)
        ss1 = np.sin(s1)
        cs1 = np.cos(s1)
        ds1 = B * ss1 * (cs1m + B / 4 * (cs1 * (- 1 + 2 * cs1m * cs1m) - B / 6 * \
            cs1m * (- 3 + 4 * ss1 * ss1) * (-3 + 4 * cs1m * cs1m)))
        s1p = s1
        s1 = dist / (b * A) + ds1

    t = su1 * ss1 - cu1 * cs1 * cb
    lat2 = np.arctan2(su1 * cs1 + cu1 * ss1 * cb, (1 - f) * np.sqrt(sa * sa + t * t))
    l2 = np.arctan2(ss1 * sb, cu1 * cs1 - su1 * ss1 * cb)
    c = f / 16 * csa * (4 + f * (4 - 3 * csa))
    l = l2 - (1 - c) * f * sa * (s1 + c * ss1 * (cs1m + c * cs1 * (-1 + 2 * cs1m * cs1m)))
    d = np.arctan2(sa, -t)
    finaltc = d + 2 * np.pi
    backtc = d + np.pi
    lon2 = lon1 + l

    return (np.rad2deg(lat2), np.rad2deg(lon2))

def shaded_great_circle(m, lat_0, lon_0, dist=100, alpha=0.2, col='k'):  # dist specified in nautical miles
    dist = dist * 1852  # Convert distance to nautical miles
    theta_arr = np.linspace(0, np.deg2rad(360), 100)
    lat_0 = np.deg2rad(lat_0)
    lon_0 = np.deg2rad(lon_0)

    coords_new = []

    for theta in theta_arr:
        coords_new.append(calc_new_coord(lat_0, lon_0, theta, dist))

    lat = [item[0] for item in coords_new]
    lon = [item[1] for item in coords_new]

    x, y = m(lon, lat)
    m.plot(x, y, col)

def get_instrument_locations(filename='WeatherInstrumentationLocation.xlsx'):
    df = pd.read_excel(filename)
    df = df[(df['IsActive'] == True) & (df['LocationType'].str.contains("Field Mill"))]
    return df[['SiteName', 'Latitude', 'Longitude', 'Elevation']].set_index('SiteName')

def make_map(ax, projection='mill', resolution='f',
             color_land='#D2B48C', color_water='#add8e6',
             llcrnrlon=-80.9, llcrnrlat=28.3 ,
             urcrnrlon=-80.45, urcrnrlat=28.8,
             ):
    m = Basemap(ax=ax, projection=projection, lon_0=-78., resolution=resolution,
                llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                )
    figure(num=None, figsize=(8, 6), dpi=300, facecolor='w', edgecolor='w')
    m.fillcontinents(color=color_land, lake_color=color_water)
    m.drawmapboundary(fill_color=color_water)
    m.drawcoastlines()
    return m

def get_readings(filename=f'field-mill-export-{day}.xlsx'):

    print(filename)
    df = pd.read_excel(filename)
    df['UTC'] = df.apply(lambda r: datetime.combine(r['Event Date'], r['Event Time']).replace(tzinfo=timezone('UTC')), 1)
    df['ET'] = df.apply(lambda r: r['UTC'].astimezone(timezone('US/Eastern')), 1)
    df.drop(['Event Date', 'Event Time'], axis=1, inplace=True)
    data = df.to_dict(orient='index')
    return df, data

def plotmap(ax, wall_clock, countdown_clock):
    m = make_map(ax)
    label_point(ax, 28.608389, -80.604333, m, label='LC39A', color='k', size=50, fontsize=10)

    # label_launch_pads(m)
    df = get_instrument_locations()

    field_mills = df.to_dict(orient='index')
    for k, v in field_mills.items():
        label_point(ax, v['Latitude'], v['Longitude'], m, marker='o')
    df_readings, data = get_readings()
    prev = None

    for k, v in data.items():
        if prev == None:
            prev = v['ET']
        fmkey = f"FM{v['Mill Number']:02d}"
        if fmkey in field_mills.keys() and v['ET'] == wall_clock:
            if abs(v['One Minute Mean']) >= threshold:
                color = 'Red'
            else:
                color = 'Green'
            label_point(ax, field_mills[fmkey]['Latitude'], field_mills[fmkey]['Longitude'],
                        m, color=color, marker='o')

    return m

def label_point(ax, lat, lon, m, color='grey', size=25, fontsize=8,
                marker='s', label=None, labeloffset=10):
    if type(lon) is not list:
        x, y = m([lon], [lat])
    else:
        x, y = m(lon, lat)
    x2, y2 = (labeloffset, 0)
    ax.scatter(x, y, size, marker=marker, color=color, zorder=2)
    if label is not None:
        ax.annotate(label, xy=(x[0], y[0]), xycoords='data',
                     xytext=(x2, y2), textcoords='offset points',
                     fontsize = fontsize, color=color)

def countdownclock(clock, t0):
    seconds = (clock - t0).total_seconds()
    min, sec = divmod(seconds, 60)
    hour, min = divmod(min, 60)
    return "T %d:%02d" % (hour, min)

class MyTestCase(unittest.TestCase):

    def test_generate_plots(self):

        df, data = get_readings()
        df = df.sort_values(by='ET', ascending=True)
        a = df['ET'].unique()
        x_values = []
        y_values = []
        oldmax = 0

        cnt = 0
        for clock in tqdm(a):
            filename = f'plots_{day}/{cnt:04}.png'
            cnt +=1
            if os.path.isfile(filename):
                continue
            minmeans = []
            for dv in data.values():
                if dv['ET'] == clock:
                    minmeans.append(abs(dv['One Minute Mean']))
            maxvalues = max(minmeans)
            oldmax = max([oldmax, maxvalues])
            x_values.append(clock)
            y_values.append(maxvalues)
            # if os.path.isfile(filename) or maxvalues < 1500:
            #     continue

            fig = plt.figure(figsize=(8, 9))
            gs = gridspec.GridSpec(2, 1, height_ratios=[8, 1])

            ax0 = plt.subplot(gs[0])
            ax1 = plt.subplot(gs[1])

            plotmap(ax0, clock, countdownclock(clock, t0))

            for spine in ["left", "top", "right", "bottom"]:
                ax1.spines[spine].set_visible(False)
            ax1.get_yaxis().set_ticks([])
            ax1.get_xaxis().set_ticks([])
            ax1.set_ylim(0, oldmax)

            ax1.axvline(x=t0)
            ax1.annotate("launch\nwindow", xytext=(0, -20), textcoords='offset points',
                         xy=(mdates.date2num(t0), 0),  ha='center')
            #
            ax1.axvline(x=limit)
            ax1.annotate("15 min\ndecision", xytext=(0, -20), textcoords='offset points',
                         xy=(mdates.date2num(limit), 0),  ha='center')

            ax1.fill_between(x_values, threshold, max(y_values), color='red', alpha=0.25)
            ax1.fill_between(x_values, 0, threshold, color='green', alpha=0.25)
            ax1.plot(x_values, y_values, color='blue')

            ax0.set_title("electric field mill launch weather criteria")
            ax1.set_title(f"{clock.strftime('%Y-%m-%d')} Demo-2 Launch Attempt")

            plt.draw()
            plt.tight_layout()
            print(x_values[-1], y_values[-1])
            fig.savefig(filename, dpi=300)
            plt.clf()



if __name__ == '__main__':
    unittest.main()
