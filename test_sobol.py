%run imports.py

sob = stats.qmc.Sobol(d=7)
gino = sob.random_base2(m=6)

stats.qmc.discrepancy(gino)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

writer = ImageMagickFileWriter(fps = 20)

with writer.saving(fig, 'sobol_3d.gif', 50):
    for angle in range(0, 360):
        ax.view_init(30, angle)
        writer.grab_frame()
