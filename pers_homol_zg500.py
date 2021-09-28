%run imports.py
cart = '/home/fedef/Research/lavori/WRtool_tests/ERA_ref/'
fi = 'out_ERA_NDJFM_EAT_4clus_4pcs_1957-2018_dtr.p'
gigi = ctl.load_wrtool(cart + fi)
gigi
gigi.keys()
plt.ion()
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111, projection='3d')
gigi['pcs'].shape
ax.scatter(*gigi['pcs'][:, :3])
ax.scatter(gigi['pcs'][:, :3])
ax.scatter(gigi['pcs'][:, 0], gigi['pcs'][:, 1], gigi['pcs'][:, 2])
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(gigi['pcs'][:, 0], gigi['pcs'][:, 1], gigi['pcs'][:, 2], s = 1)
fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot(111, projection='3d')
for reg in range(4):
    okreg = gigi['labels'] == reg
    ax.scatter(gigi['pcs'][okreg, 0], gigi['pcs'][okreg, 1], gigi['pcs'][okreg, 2], s = 1)
import gudhi
gudhi.RipsComplex?
pino = gudhi.RipsComplex(gigi['pcs'])
pino.create_simplex_tree?
prot = pino.create_simplex_tree()
prot
prot.compute_persistence()
prot.write_persistence_diagram()
prot.write_persistence_diagram?
prot.write_persistence_diagram(open('prova_pers.txt', 'wb'))
prot.write_persistence_diagram('prova_pers.txt')
