import sys
import os
import numpy as np
import json


import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
sns.set_style('white')
sns.set_context('paper')

js = [json.load(open(f)) for f in sys.argv[1:]]

methods = [os.path.dirname(f).split("/")[1] for f in sys.argv[1:]]

df = pd.DataFrame.from_records(js)
df['method'] = methods
df['FDR'] = 1 - df['precision']
df.to_csv("supp-table1.tsv", sep="\t", index=False, float_format="%.3f")
print("wrote supp-table1.tsv")

sdf = df['method FDR FN   FP  TP-call precision recall f1'.split()]
m = pd.melt(sdf, id_vars=['method'], value_vars=['recall', 'FDR'])


#ax = sns.barplot(x="method", hue="variable", y="value", data=m)
#ax.set_ylim(0.0, 0.9)

fig , ax = plt.subplots(1, 1, sharex=True)

#ax[0].bar(1 + np.arange(len(js)), df['FP'])
#ax[1].bar(1 + np.arange(len(js)), df['TP-call'])

#ax[1].set_xticks(1 + np.arange(len(js)))
#ax[1].set_xticklabels(methods)

ax.scatter(df['FP'], df['TP-call'])
ax.set_xlim(0)
ax.set_ylim(0)
ax.set_xlabel("False Positives")
ax.set_ylabel("True Positives")

lookup = {"dhbfc": "DHBFC < 0.75", "dhd": "DHD < 0",
         "both": "(DHBFC < 0.75) && (DHD < 0)",
         "stringent": "(DHBFC < 0.65) && (DHD == -2)",
         "unfiltered": "unfiltered"}

pos = ['baseline', 'top']
hpos = ['center', 'left', 'right', 'right', 'center']
for i, txt in enumerate(df['method']):
    ax.annotate(lookup[txt], (df['FP'][i], df['TP-call'][i]),
            horizontalalignment=hpos[i], verticalalignment=pos[i%2],
            fontsize=11)

sns.despine()
plt.show()
