latbins = np.arange(-90, 91, 20)
lacol = ctl.color_set(9, use_seaborn=False)
for ru, ax in zip(allru2, axs.flatten()):
    coso = yeamean[(ru, var)]
    coso1 = coso[:20].mean(axis = 0)
    cosoanom = coso-coso1
    ax.set_title(ru)

    for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
        print(la1,la2)
        cosolat = cosoanom.sel(lat = slice(la1, la2)).mean(['lat','lon'])
        cosmu = ctl.butter_filter(cosolat, 50)
        ax.plot(coso.year, cosmu, color = col)

fig, axs = plt.subplots(2, 3, figsize = (16,9))
for ru, ax in zip(allru2, axs.flatten()):
    coso = yeamean[(ru, var)]
    coso1 = coso[:20].mean(axis = 0)
    cosoanom = coso-coso1
    ax.set_title(ru)
    smut = 50
    if ru in ['ssp585', 'hist']:
        smut = 20

    for la1, la2, col in zip(latbins[:-1], latbins[1:], lacol):
        print(la1,la2)
        cosolat = cosoanom.sel(lat = slice(la1, la2)).mean(['lat','lon'])
        cosmu = ctl.butter_filter(cosolat, smut)
        ax.plot(coso.year, cosmu, color = col)
