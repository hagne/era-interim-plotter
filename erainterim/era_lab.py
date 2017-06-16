import xarray as _xr
import matplotlib.pylab as _plt
from mpl_toolkits.basemap import Basemap as _Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable
import matplotlib.colors as _colors
import numpy as _np
import pandas as _pd

# era_interim_params = [{'name_internal':'wind_u10',
#                        'name_nc_variable': 'u10'},
#                       {'name_internal':'wind_v10',
#                        'name_nc_variable': 'v10'},
#                       {'name_internal':'wind_speed10',
#                        'name_nc_variable': None}]


def load_nc(fname):
    #     nc = Dataset(fname)
    xr_set = _xr.open_dataset(fname)
    era_int = EraInterim()
    for var in xr_set.variables.keys():
        if var in ['longitude', 'latitude', 'time']:
            setattr(era_int, var, xr_set.variables[var])
        else:
            setattr(era_int, "_" + var, xr_set.variables[var])

    era_int._generate_properties()
    return era_int

class EraInterim(object):
    def _generate_properties(self):
        if hasattr(self, '_v10') and hasattr(self, '_u10'):
            self.wind_10 = Wind(self, self._v10, self._u10)


class AtmProperty(object):
    def plot_map(self, extend=1000000, lat_0=42.030484, lon_0=-70.049392, ax = None):
        if ax:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()
        bmap = _Basemap(projection='aeqd',
                        lat_0=lat_0,
                        lon_0=lon_0,
                        width=extend,
                        height=extend,
                        resolution='i',
                        ax=a)
        # bmap.shadedrelief()
        bmap.bluemarble()

        bmap.drawcoastlines(zorder=100)
        return f, a, bmap

class Wind(AtmProperty):
    def __init__(self, era_int, v, u):
        self._era_int = era_int
        self.v = v
        self.u = u

        self._wind_speed = None

    @property
    def wind_speed(self):
        if type(self._wind_speed).__name__ == 'NoneType':
            self._wind_speed = _np.sqrt(self.u ** 2 + self.v ** 2)
        return self._wind_speed


    def plot_streamlines(self, time = 'all', vmin = None, vmax = None, ax = None):
        """
        Parameters
        ----------
        time: string, int, or Timestamp (Timestamp not working right now)
            When more then one dataset is included you can select a particular by eigher giving the the index (int) or
            the timestamp (pandas Timestamp or string). Or you can choose to plot all in separate subplots.
        """
        if type(time) == int:
            f, a, sl, cb = self._plot_streamlines(time, vmin = vmin, vmax = vmax, ax = ax)
        elif time == 'all':
            if type(ax) != type(None):
                raise ValueError('Ax has to be None if multiple plots are beeing plotted. Select a particular plot or set ax to None.')
            f,a = _plt.subplots(1,len(self._era_int.time))
            for e, at in enumerate(a):
                f, at, sl, cb = self._plot_streamlines(e, ax = at, vmin = vmin, vmax = vmax)

        return f,a,sl,cb


    def _plot_streamlines(self, time2use_idx, ax = None, vmin = None, vmax = None):
        f, a, bmap = super().plot_map(ax = ax)

        # generate evenly spaced grid
        nop = 40
        lons, lats, x, y = bmap.makegrid(nop, nop, returnxy=True)

        # create meshgrid from lon and lat
        v10psel = self.v[time2use_idx, :, :]
        u10psel = self.u[time2use_idx, :, :]

        X, Y = _np.meshgrid(self._era_int.longitude, self._era_int.latitude)

        # apply projection
        xm_lon, ym_lat = bmap(X, Y)

        # project to evenly speced grid
        new_v10psel = _plt.mlab.griddata(xm_lon.flatten(), ym_lat.flatten(), v10psel.values.flatten(), x, y,
                                         interp='linear')
        new_u10psel = _plt.mlab.griddata(xm_lon.flatten(), ym_lat.flatten(), u10psel.values.flatten(), x, y,
                                         interp='linear')
        speed = _np.sqrt(new_u10psel ** 2 + new_v10psel ** 2)

        # plot streamlines
        cmap = _plt.cm.YlOrBr

        lwmax = 3
        lwmin = 0.2
        lw = (speed * (lwmax - lwmin) / (vmax - vmin))
        lw += lwmin - lw.min()

        sl = a.streamplot(x, y, new_v10psel, new_u10psel,
                          #                   linewidth=1,
                          density=3,
                          color=speed,
                          cmap=cmap,
                          norm=_colors.Normalize(vmin=vmin, vmax=vmax),
                          linewidth=lw  # 0.5*speed
                          )
        # a.pcolormesh(x,y,new_z)
        divider = _make_axes_locatable(a)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb = f.colorbar(sl.lines, cax=cax)
        cb.set_label('Wind speed m/s')

        time2use = _pd.Timestamp(self._era_int.time.values[time2use_idx])
        a.set_title(time2use, loc='left')
        return f, a, sl, cb
