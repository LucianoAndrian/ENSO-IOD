#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer

#esto deberia poder hacerse en una linea
months = []
for num in range(1900, 2011, 1):
  months.append(('/'.join([f'{num}' '%.2d' '01' % m for m in range(1,13)])))

months = ('/'.join(months))

server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="ba9ad8fcbb4f9f895b228f546713d1ba",email="luciano.andrian@cima.fcen.uba.ar")
server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": months,
    "expver": "1",
    "levtype": "sfc",
    "param": "151.128",
    "stream": "moda",
    "type": "an",
    "grid": "0.25/0.25",
    "format": "netcdf",
    "target": "ERA-20C.mon.psl.nc",
})

server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="ba9ad8fcbb4f9f895b228f546713d1ba",email="luciano.andrian@cima.fcen.uba.ar")
server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": months,
    "expver": "1",
    "levtype": "PL",
    "levelist":"200/500",
    "param": "129.128",
    "stream": "moda",
    "domain": "G",
    "type": "an",
    "grid": "0.25/0.25",
    "format": "netcdf",
    "target": "ERA-20C.geopot500_200.mon.nc",
})

server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="ba9ad8fcbb4f9f895b228f546713d1ba",email="luciano.andrian@cima.fcen.uba.ar")
server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": months,
    "expver": "1",
    "levtype": "PL",
    "levelist":"500",
    "param": "135.128",
    "stream": "moda",
    "domain": "G",
    "type": "an",
    "grid": "0.25/0.25",
    "format": "netcdf",
    "target": "ERA-20C.mon.w.mon.nc",
})


server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="ba9ad8fcbb4f9f895b228f546713d1ba",email="luciano.andrian@cima.fcen.uba.ar")
server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": months,
    "expver": "1",
    "levtype": "PL",
    "levelist":"200",
    "param": "155.128",
    "stream": "moda",
    "domain": "G",
    "type": "an",
    "grid": "0.25/0.25",
    "format": "netcdf",
    "target": "ERA-20C.mon.div200.nc",
})

server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="ba9ad8fcbb4f9f895b228f546713d1ba",email="luciano.andrian@cima.fcen.uba.ar")
server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": months,
    "expver": "1",
    "levtype": "PL",
    "levelist":"200",
    "param": "131.128/132.128",
    "stream": "moda",
    "domain": "G",
    "type": "an",
    "grid": "0.25/0.25",
    "format": "netcdf",
    "target": 'ERA-20C.mon.U-V_200_50-78.nc',
})


server = ECMWFDataServer(url="https://api.ecmwf.int/v1",key="ba9ad8fcbb4f9f895b228f546713d1ba",email="luciano.andrian@cima.fcen.uba.ar")
server.retrieve({
    "class": "e2",
    "dataset": "era20c",
    "date": months,
    "expver": "1",
    "levtype": "sfc",
    "param": "167.128",
    "stream": "moda",
    "domain": "G",
    "type": "an",
    "grid": "0.5/0.5",
    "format": "netcdf",
    "target": 'pp_ERA20c.nc',
})












