#!/usr/bin/env python

import sys
import glob
import pandas as pd
import pyarrow.parquet as pq


if len(sys.argv)>2:
	chrc = sys.argv[2].split(":")[0]
	startc = int(sys.argv[2].split(":")[1].split("-")[0])
	endc = int(sys.argv[2].split(":")[1].split("-")[1])
	#print(endc)
	pq.read_table(sys.argv[1],filters=[('chrom','=',chrc),('pos','>=',startc),('pos','<=',endc)]).to_pandas().to_csv(sys.stdout,index=False,sep="\t")
else:
	pd.read_parquet(sys.argv[1]).to_csv(sys.stdout,index=False,sep="\t")

