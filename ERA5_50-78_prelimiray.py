import cdsapi

c = cdsapi.Client()

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
    },
    'ERA5_slp_50-78.nc')

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'variable': 'geopotential',
        'pressure_level': [
            '200',
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
    },
    'ERA5_HGT200_50-78.nc')


c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'variable': 'vertical_velocity',
        'pressure_level': '500',
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
    },
    'ERA5_vp500_50-78.nc')


c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'variable': 'vertical_velocity',
        'pressure_level': '200',
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
    },
    'ERA5_vp200_50-78.nc')


c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'variable': 'divergence',
        'product_type': 'reanalysis-monthly-means-of-daily-means',
        'pressure_level': '200',
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
    'ERA5.mon.div200_50-78.nc')


import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means-preliminary-back-extension',
    {
        'format': 'netcdf',
        'product_type': 'reanalysis-monthly-means-of-daily-means',
        'variable': [
            'u_component_of_wind', 'v_component_of_wind',
        ],
        'pressure_level': '200',
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
    'ERA5.mon.U-V_200_50-78.nc')