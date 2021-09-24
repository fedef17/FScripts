
%run bottino_overview_5.py

glomeans
plt.ion()
coso = yeamean[(ru, var)]
ru = 'b050'
var = 'tas'
coso = yeamean[(ru, var)]
coso
coso = glomeans[(ru, var)]
coso
coso.plot()
plt.plot(*coso)
coef2, covmat2 = np.polyfit(coso[0], coso[1], deg = 2, cov = True)
coef3, covmat3 = np.polyfit(coso[0], coso[1], deg = 3, cov = True)
years = coso[0]
fitco2 = np.polyval(coef2, years)
fitco3 = np.polyval(coef3, years)
plt.plot(years, fitco2)
plt.plot(years, fitco3)
plt.figure()
plt.plot(years, coso-fitco2)
plt.plot(years, coso[1]-fitco2)
plt.plot(years, coso[1]-fitco3)
plt.figure()
ps = np.abs(np.fft.rfft(coso[1]))**2
frq = np.fft.rfftfreq(coso[1].size, 1)
plt.plot(frq, ps)
gtas = coso[1]
gtas
gtas = coso[1] - np.mean(coso[1])
gtas
ps = np.abs(np.fft.rfft(gtas))**2
plt.figure()
plt.plot(frq, ps)
ps2 = np.abs(np.fft.rfft(fitco2))**2
ps3 = np.abs(np.fft.rfft(fitco3))**2
plt.plot(frq, ps2)
plt.plot(frq, ps3)
ps2 = np.abs(np.fft.rfft(fitco2 - np.mean(fitco2)))**2
ps3 = np.abs(np.fft.rfft(fitco3 - np.mean(fitco3)))**2
plt.figure()
plt.plot(frq, ps)
plt.plot(frq, ps2)
plt.plot(frq, ps3)
invfr = 1/frq
plt.figure()
plt.plot(invfr, ps)
plt.plot(invfr, ps2)
plt.plot(invfr, ps3)
ps2 = np.abs(np.fft.rfft(fitco2))**2
fitco2
ps2 = np.abs(np.fft.rfft(coso[1] - fitco2))**2
ps3 = np.abs(np.fft.rfft(coso[1] - fitco3))**2
plt.figure()
plt.plot(invfr, ps)
plt.plot(invfr, ps2)
plt.plot(invfr, ps3)
gtas2 = coso[1] - fitco2
gtas3 = coso[1] - fitco3
tama = yeamean[(ru, var)]
fitco2
tama2 = tama - fitco2
tama2 = tama - fitco2[:, np.newaxis, np.newaxis]
tama3 = tama - fitco3[:, np.newaxis, np.newaxis]
gtas2
g50_2 = ctl.butter_filter(gtas2, 50)
g50_3 = ctl.butter_filter(gtas3, 50)
plt.figure(2)
plt.plot(years, g50_2, color = 'orange', linewidth = 3)
plt.plot(years, g50_3, color = 'orange', linewidth = 3)
plt.plot(years, g50_2, color = 'blue', linewidth = 3)
plt.figure()
plt.plot(years, np.gradient(g50_2))
plt.plot(years, np.gradient(g50_3))
plt.grid()
incr = (years > 2100) & (np.gradient(g50_3) > 0)
decr = (years > 2100) & (np.gradient(g50_3) < 0)
tasincr = tama2[incr]
tasdecr = tama3[decr]
tasincr = tama3[incr]
ctl.plot_multimap_contour([tasincr, tasdecr])
ctl.plot_map_contour(tasincr)
tasincr = tama3[incr].mean('year')
tasdecr = tama3[decr].mean('year')
ctl.plot_multimap_contour([tasincr, tasdecr])
tasdecr = tama3[decr].mean('year') - tama3.mean('year')
tasincr = tama3[incr].mean('year') - tama3.mean('year')
ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9))
ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9), plot_anomalies=True)
ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9), plot_anomalies=True, subtitles= ['gtas increasing', 'gtas decreasing'])
ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9), plot_anomalies=True, subtitles= ['gtas increasing', 'gtas decreasing'], cbar_range=(-0.5, 0.5))
tasincr = tama3[incr].mean('year') - tama3[years > 2100].mean('year')
tasdecr = tama3[decr].mean('year') - tama3[years > 2100].mean('year')
ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9), plot_anomalies=True, subtitles= ['gtas increasing', 'gtas decreasing'], cbar_range=(-0.5, 0.5))
ctl.plot_multimap_contour([tasincr, tasdecr], figsize = (16,9), plot_anomalies=True, subtitles= ['gtas increasing', 'gtas decreasing'], cbar_range=(-0.3, 0.3))
