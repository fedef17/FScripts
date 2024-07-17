# coding: utf-8
get_ipython().run_line_magic('run', 'imports.py')
cart = '/home/fedef/Research/lavori/WRtool_tests/ERA_ref/'
fi = 'out_ERA_NDJFM_EAT_4clus_4pcs_1957-2018_dtr.p'
gigi = ctl.load_wrtool(cart + fi)
gigi
gigi['pcs']
gigi['pcs'].shape
9211/90
get_ipython().run_line_magic('pinfo', 'ctl.seasonal_set')
ctl.seasonal_set(gigi['pcs'], dates = gigi['dates'], season = 'NDJFM')
zupi = ctl.seasonal_set(gigi['pcs'], dates = gigi['dates'], season = 'NDJFM')
zupi
len(zupi)
zupik, dates_zupik = zupi
zupik
zupik.shape
zupik_low = []
for zup in zupik:
    zupik_low.append(ctl.running_mean(zup, 5))
    
zupik_low = np.stack(zupik_low)
zupik_low.shape
reload(ctl)
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik)
reload(cd)
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik)
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik, flag_perc=True)
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik, flag_perc=True, numclus = 5)
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik, flag_perc=True, perc = 70, numclus = 5)
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik, flag_perc=True, perc = 70, numclus = 4)
centroids
centroids.shape
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik, numpcs=3, numclus = 4)
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
labels
colors = ctl.color_set(4)
plt.ion()
plt.show()
centroids, labels, cluspatterns, repres, distances = cd.EnsClus_light(zupik_low, numpcs=3, numclus = 4)
for zup, lab in zip(zupik_low, labels):
    ax.scatter(zup[:,0], zup[:,1], zup[:,2], color = colors[lab], s = 1)
    ax.plot(zup[:,0], zup[:,1], zup[:,2], color = colors[lab], alpha = 0.5)
    
ax.clear()
for zup, lab in zip(zupik_low, labels):
    #ax.scatter(zup[:,0], zup[:,1], zup[:,2], color = colors[lab], s = 1)
    ax.plot(zup[:,0], zup[:,1], zup[:,2], color = colors[lab], alpha = 0.5)
    
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
for nu, col in zip([17,48,24,8], colors):
    pio = zupik_low[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], color = col, alpha = 0.7)
    
zupik_diff = np.diff(zupik_low, axis = 1)
zupik_diff.shape
centroids_diff, labels_diff, cluspatterns_diff, repres_diff, distances_diff = cd.EnsClus_light(zupik_diff, numpcs=3, numclus = 4)
centroids_diff, labels_diff, cluspatterns_diff, repres_diff, distances_diff = cd.EnsClus_light(zupik_diff, numpcs=10, numclus = 4)
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
for nu, col in zip(np.arange(4), colors):
    pio = cluspatterns[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], color = col, alpha = 0.7)
    
    
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
for nu, col in zip(np.arange(4), colors):
    pio = cluspatterns_diff[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], color = col, alpha = 0.7)
    
    
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
for nu, col in zip(repres_diff, colors):
    pio = zupik_low[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], color = col, alpha = 0.7)
    
    
repres_diff
repres
for nu, col in zip(repres_diff, colors):
    pio = zupik_low[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], color = col, alpha = 0.7)
    
    
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
zupik_diff
zupik_diff.shape
zupik_speed = np.sum(zupik_diff**2, axis = -1)
zupik_speed
zupik_speed.shape
centroids_speed, labels_speed, cluspatterns_speed, repres_speed, distances_speed = cd.EnsClus_light(zupik_speed, numpcs=10, numclus = 4)
centroids_speed, labels_speed, cluspatterns_speed, repres_speed, distances_speed = cd.EnsClus_light(zupik_speed, numpcs=4, numclus = 4)
for nu, col in zip(repres_speed, colors):
    pio = zupik_low[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], color = col, alpha = 0.7)
    
    
    
fig = plt.figure(figsize=(8, 6), dpi=150)
ax = fig.add_subplot(111, projection='3d')
for nu, col in zip(repres_speed, colors):
    pio = zupik_low[nu]
    ax.plot(pio[:, 0], pio[:,1], pio[:,2], alpha = 0.7)
    
    
    
    
centroids_speed, labels_speed, cluspatterns_speed, repres_speed, distances_speed = cd.EnsClus_light(zupik_speed, numpcs=4, numclus = 3)
centroids_speed, labels_speed, cluspatterns_speed, repres_speed, distances_speed = cd.EnsClus_light(zupik_speed, numpcs=4, numclus = 4)
centroids_speed, labels_speed, cluspatterns_speed, repres_speed, distances_speed = cd.EnsClus_light(zupik_speed, numpcs=10, numclus = 4)
centroids_speed
plt.figure()
cluspatterns_speed
cluspatterns_speed.shape
for co in cluspatterns_speed:
    plt.plot(co)
    
