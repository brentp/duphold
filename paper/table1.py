import sys
import os
import numpy as np
import json


import pandas as pd
import os
import sys
import json
import tabulate

js = [json.load(open(f)) for f in sys.argv[1:]]

methods = [os.path.dirname(f).split("/")[1] for f in sys.argv[1:]]

df = pd.DataFrame.from_records(js)
df['method'] = methods
df['FDR'] = 1 - df['precision']
df.to_csv("supp-table1.tsv", sep="\t", index=False, float_format="%.3f")
print("wrote supp-table1.tsv")

sdf = df['method FDR FN   FP  TP-call precision recall f1'.split()]
sdf.to_csv("table1.tsv", sep="\t", index=False, float_format="%.3f")
print("wrote table1.tsv")

sdf = sdf.drop('f1', axis=1)
sdf['recall-%'] = 100.0 * sdf['recall'] / float(sdf['recall'][0])
sdf['FP-%'] = 100.0 * sdf['FP'] / float(sdf['FP'][0])
print(tabulate.tabulate(sdf.values, sdf.columns, tablefmt="pipe", floatfmt=".3f"))

