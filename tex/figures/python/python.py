import matplotlib.pyplot as pl
fig, ax = pl.subplots(1)
ax.annotate("Placeholder", xy=(0.5, 0.5), xytext=(0, 0),
            xycoords='axes fraction', textcoords='offset points',
            ha='center', va='center', fontsize=20)
ax.axis('off')
fig.savefig("python.pdf", bbox_inches='tight')
