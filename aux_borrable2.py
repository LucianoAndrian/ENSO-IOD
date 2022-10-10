sst = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
sst = sst.sel(time=slice('1970-01-01', '2020-12-31'))

coslat = np.cos(np.deg2rad(sst.coords['lat'].values))
wgts = np.sqrt(coslat)[..., np.newaxis]

solver = Eof(sst['sst'], weights=wgts, center=False)


eof_sst = solver.eofs(neofs=5)
cps_sst = solver.pcs(npcs=5)

var1_obs, var2_obs, var3_obs, var4_obs, var5_obs = np.around(solver.varianceFraction(neigs=5).values*100,1)
eoc_cor_obs = solver.eofsAsCorrelation(neofs=3)

plt.plot(cps_sst.sel(mode=1));plt.show()
plt.imshow(eof_sst.sel(mode=1),cmap='RdBu',vmax=0.05, vmin=-0.05);plt.colorbar();plt.show()