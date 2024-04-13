import urllib.request
from bs4 import BeautifulSoup
import netCDF4
import argparse
import datetime
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature
from cartopy.util import add_cyclic_point
import matplotlib.path as mpath



def descarga_nc(inid, inim, find, finm, finy, var_desc,var, level):
    # Open NCEP NCAR access to link to data
    url = 'http://www.psl.noaa.gov/cgi-bin/data/composites/comp.day.pl?var='+var_desc+'&level='+level+'&iy[1]=&im[1]=&id[1]=&iy[2]=&im[2]=&id[2]=&iy[3]=&im[3]=&id[3]=&iy[4]=&im[4]=&id[4]=&iy[5]=&im[5]=&id[5]=&iy[6]=&im[6]=&id[6]=&iy[7]=&im[7]=&id[7]=&iy[8]=&im[8]=&id[8]=&iy[9]=&im[9]=&id[9]=&iy[10]=&im[10]=&id[10]=&iy[11]=&im[11]=&id[11]=&iy[12]=&im[12]=&id[12]=&iy[13]=&im[13]=&id[13]=&iy[14]=&im[14]=&id[14]=&iy[15]=&im[15]=&id[15]=&iy[16]=&im[16]=&id[16]=&iy[17]=&im[17]=&id[17]=&iy[18]=&im[18]=&id[18]=&iy[19]=&im[19]=&id[19]=&iy[20]=&im[20]=&id[20]&monr1='+str(inim)+'&dayr1='+str(inid)+'&monr2='+str(finm)+'&dayr2='+str(find)+'&iyr[1]='+str(finy)+'&filenamein=&plotlabel=&lag=0&labelc=Color&labels=Shaded&type=2&scale=&label=0&cint=&lowr=&highr=&istate=0&proj=ALL&xlat1=&xlat2=&xlon1=&xlon2=&custproj=Cylindrical+Equidistant&level1=1000mb&level2=10mb&Submit=Create+Plot'
    response = urllib.request.urlopen(url)
    data = response.read()      # a `bytes` object
    soup = BeautifulSoup(data,'html.parser') #is an xml, beautifull has a module to manage it
	#A very inefficient way to get the nc file
    link = soup.findAll('img')[-1]['src']
    link=list(link)
    link[-3]='n'
    link[-2]='c'
    link[-1]=''
    #get nc file save as netcdf
    ruta = "./tmp/"
    url_nc='http://www.psl.noaa.gov'+"".join(link)
    urllib.request.urlretrieve(url_nc, ruta+var+'.nc')

descarga_nc(inid, inim, find, finm, finy, var_desc,var, level)
