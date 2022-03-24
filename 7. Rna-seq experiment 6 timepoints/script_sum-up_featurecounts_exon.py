#!/usr/bin/env python
import pandas as pd
counts = pd.read_table('counts.Tfas_Tlei_6_timepoints.exons.edited.forR.txt')
print(counts.head())
counts = counts.groupby('Geneid').sum()
counts.to_csv(r'counts.Tfas_Tlei_6_timepoints.exons.edited.forR.sum.txt', index=True, sep='\t', mode='a')
