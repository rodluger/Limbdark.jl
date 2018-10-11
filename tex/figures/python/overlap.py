import matplotlib.pyplot as pl
import numpy as np
from matplotlib.patches import Arc

edgecolor = 'k'
emittercolor = 'k'
occultorcolor = 'k'
facecolor = '#7cbce9'
r = 0.75
b = 1.25

phi = np.linspace(0, 2 * np.pi, 1000)
x = np.cos(phi)
y = np.sin(phi)

xo = r * np.cos(phi)
yo = b + r * np.sin(phi)

fig, ax = pl.subplots(1, 4, figsize=(12, 5))
for axis in ax:
    axis.plot(x, y, '-', color=emittercolor, lw=1, zorder=0)
    axis.plot(xo, yo, '-', color=occultorcolor, lw=1, zorder=0)
    axis.set_xlim(-1.25, 1.25)
    axis.set_ylim(-1.25, 2.25)
    axis.set_aspect(1)
    axis.axis('off')

x = np.linspace(-0.6, 0.6, 1000)
y1 = np.sqrt(1 - x ** 2)
y2 = b - np.sqrt(r ** 2 - x ** 2)
ax[0].fill_between(x, y1, y2, zorder=1, edgecolor=edgecolor, facecolor=facecolor)
ax[0].set_title("area of overlap")

y1 = np.concatenate((np.linspace(0.8, b, 500), np.linspace(b, 0.8, 500)))
y2 = b - np.sqrt(r ** 2 - x ** 2)
ax[1].fill_between(x, y1, y2, zorder=1, edgecolor=edgecolor, facecolor=facecolor)
ax[1].set_title("occultor sector")

y1 = np.sqrt(1 - x ** 2)
y2 = np.concatenate((np.linspace(0.8, 0, 500), np.linspace(0, 0.8, 500)))
ax[2].fill_between(x, y1, y2, zorder=1, edgecolor=edgecolor, facecolor=facecolor)
ax[2].set_title("emitter sector")

y1 = np.concatenate((np.linspace(0.8, b, 500), np.linspace(b, 0.8, 500)))
y2 = np.concatenate((np.linspace(0.8, 0, 500), np.linspace(0, 0.8, 500)))
ax[3].fill_between(x, y1, y2, zorder=1, edgecolor=edgecolor, facecolor=facecolor)
ax[3].set_title("kite area")

for y in [0.5, 0.785]:
    pl.figtext(0.3, y, "=", fontsize=18)
    pl.figtext(0.5, y, "+", fontsize=16)
    pl.figtext(0.71, y, "-", fontsize=20)


ax[1].plot([0, 0], [0, b], 'k-', lw=1, zorder=2)
ax[1].plot([0, 0.6], [0, 0.8], 'k-', lw=1, zorder=2)
ax[1].annotate(r"$r$", xy=(0.34, 1.13), xycoords="data", xytext=(0, 0),
               textcoords="offset points", ha="center", va="center",
               fontsize=14)
ax[1].annotate(r"$1$", xy=(0.34, 0.29), xycoords="data", xytext=(0, 0),
               textcoords="offset points", ha="center", va="center",
               fontsize=14)
ax[1].annotate(r"$b$", xy=(-0.12, 0.69), xycoords="data", xytext=(0, 0),
               textcoords="offset points", ha="center", va="center",
               fontsize=14)

arc1 = Arc((0, 0), 0.5, 0.5, 0, 53.3, 90)
arc2 = Arc((0, b), 0.5, 0.5, 270, 0, 53.3)
ax[1].add_artist(arc1)
ax[1].add_artist(arc2)

ax[1].annotate(r"$\kappa_0$", xy=(0.14, 0.90), xycoords="data", xytext=(0, 0),
               textcoords="offset points", ha="center", va="center",
               fontsize=14)
ax[1].annotate(r"$\kappa_1$", xy=(0.11, 0.38), xycoords="data", xytext=(0, 0),
               textcoords="offset points", ha="center", va="center",
               fontsize=14)

fig.savefig("overlap.pdf", bbox_inches='tight')
