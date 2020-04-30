import pandas as pd
from matplotlib import rcParams
from matplotlib import style
import matplotlib.pyplot as plt

style.use('fivethirtyeight')
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Trebuchet MS']
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['xtick.color'] = 'grey'
plt.rcParams['ytick.color'] = 'grey'

df = pd.read_csv('allele_len.tsv', sep = '\t')

fig, ax = plt.subplots(figsize=(10,3))
df['LENGTH'].hist(bins = 50)
ax.set_yscale('log')
ax.set_title('TCR/BCR allele length distribution')
ax.set_xlabel('bases')
ax.set_ylabel('counts')
plot.show()
plt.show()
