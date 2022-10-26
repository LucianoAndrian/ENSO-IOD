# Descarga para TP AA.
# ERA%
# # U 850 - 200 1950-1979
# hgt 850 - 200 1950-1979,1979-2020
# psl 1950-1978,1979-2020

########################################################################################################################
import cdsapi
c = cdsapi.Client()

# U 850 - 200 1950-1979
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'product_type': 'reanalysis-monthly-means-of-daily-means',
        'variable': [
            'u_component_of_wind', 'v_component_of_wind',
        ],
        'pressure_level': ['200'],
        'year': [
            '1950', '1951', '1952',
            '1953', '1954', '1955',
            '1956', '1957', '1958',
            '1959', '1960', '1961',
            '1962', '1963', '1964',
            '1965', '1966', '1967',
            '1968', '1969', '1970',
            '1971', '1972', '1973',
            '1974', '1975', '1976',
            '1977', '1978',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',

    },
    'ERA5.mon.UV200_50-78.nc')


# hgt 850 - 200 1950-1978
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'variable': 'geopotential',
        'pressure_level': [
            '200', '850',
        ],
        'year': [
            '1950', '1951', '1952',
            '1953', '1954', '1955',
            '1956', '1957', '1958',
            '1959', '1960', '1961',
            '1962', '1963', '1964',
            '1965', '1966', '1967',
            '1968', '1969', '1970',
            '1971', '1972', '1973',
            '1974', '1975', '1976',
            '1977', '1978',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'product_type': 'reanalysis-monthly-means-of-daily-means',
        "grid": "2/2"
    },
    'ERA5.geopo850_200.mon_50-78.nc')

# hgt 1979-2020
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'geopotential',
        'pressure_level': [
            '200', '850'
        ],
        'year': [
            '1979', '1980', '1981',
            '1982', '1983', '1984',
            '1985', '1986', '1987',
            '1988', '1989', '1990',
            '1991', '1992', '1993',
            '1994', '1995', '1996',
            '1997', '1998', '1999',
            '2000', '2001', '2002',
            '2003', '2004', '2005',
            '2006', '2007', '2008',
            '2009', '2010', '2011',
            '2012', '2013', '2014',
            '2015', '2016', '2017',
            '2018', '2019', '2020',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        "grid": "2/2"
    },
    'ERA5.geopot850_200.mon.nc')


# psl 1950-1978
c.retrieve(
    'reanalysis-era5-single-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'product_type': 'reanalysis-monthly-means-of-daily-means',
        'variable': 'mean_sea_level_pressure',
        'year': [
            '1950', '1951', '1952',
            '1953', '1954', '1955',
            '1956', '1957', '1958',
            '1959', '1960', '1961',
            '1962', '1963', '1964',
            '1965', '1966', '1967',
            '1968', '1969', '1970',
            '1971', '1972', '1973',
            '1974', '1975', '1976',
            '1977', '1978',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        "grid": "2/2"
    },
    'ERA5.mon.slp_50-78.nc')

# psl 1979-2020
c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'mean_sea_level_pressure',
        'year': [
            '1979', '1980', '1981',
            '1982', '1983', '1984',
            '1985', '1986', '1987',
            '1988', '1989', '1990',
            '1991', '1992', '1993',
            '1994', '1995', '1996',
            '1997', '1998', '1999',
            '2000', '2001', '2002',
            '2003', '2004', '2005',
            '2006', '2007', '2008',
            '2009', '2010', '2011',
            '2012', '2013', '2014',
            '2015', '2016', '2017',
            '2018', '2019', '2020',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        "grid": "2/2"
    },
    'ERA5.mon.slp.nc')

########################################################################################################################

# hgt 2021
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'geopotential',
        'pressure_level': [
            '200', '850'
        ],
        'year': [
            '2021',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        "grid": "2/2"
    },
    'ERA5.geopot850_200_2021.mon.nc')



# psl 1979-2020
c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': 'mean_sea_level_pressure',
        'year': [
            '2021',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        "grid": "2/2"
    },
    'ERA5.mon.slp_2021.nc')

# U 850 - 200 1950-1979
c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            'u_component_of_wind', 'v_component_of_wind'
        ],
        'pressure_level': ['200'],
        'year': [
            '2021',
        ],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
    },
    'ERA5.mon.U200-850_2021.nc')