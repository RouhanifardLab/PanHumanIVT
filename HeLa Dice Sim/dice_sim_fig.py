import os
import pandas as pd
import numpy as np
import seaborn as sns

def d_sim(intersection, n1, n2):
    return 2 * intersection / (n1 + n2)

df = pd.read_csv('trip_intersect.csv')
df.set_index('Unnamed: 0', inplace=True)
tripsect_overlap_dict_out = df.to_dict('index')

ivt1 = {'30': (tripsect_overlap_dict_out['30%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['30%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['30%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['30%']['ivt_1']) ,
        '40': (tripsect_overlap_dict_out['40%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['40%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['40%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['40%']['ivt_1']),
        '50': (tripsect_overlap_dict_out['50%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['50%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['50%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['50%']['ivt_1']), 
        '60': (tripsect_overlap_dict_out['60%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['60%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['60%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['60%']['ivt_1']),
        '70': (tripsect_overlap_dict_out['70%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['70%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['70%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['70%']['ivt_1']), 
        '80': (tripsect_overlap_dict_out['80%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['80%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['80%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['80%']['ivt_1']),
        '95': (tripsect_overlap_dict_out['95%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['95%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['95%']['ivt_1_3'] + 
               tripsect_overlap_dict_out['95%']['ivt_1'])}
ivt2 = {'30': (tripsect_overlap_dict_out['30%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['30%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['30%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['30%']['ivt_2']),
        '40': (tripsect_overlap_dict_out['40%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['40%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['40%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['40%']['ivt_2']),
        '50': (tripsect_overlap_dict_out['50%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['50%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['50%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['50%']['ivt_2']), 
        '60': (tripsect_overlap_dict_out['60%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['60%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['60%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['60%']['ivt_2']),
        '70': (tripsect_overlap_dict_out['70%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['70%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['70%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['70%']['ivt_2']), 
        '80': (tripsect_overlap_dict_out['80%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['80%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['80%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['80%']['ivt_2']),
        '95': (tripsect_overlap_dict_out['95%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['95%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['95%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['95%']['ivt_2'])}
ivt3 = {'30': (tripsect_overlap_dict_out['30%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['30%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['30%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['30%']['ivt_2']),
        '40': (tripsect_overlap_dict_out['40%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['40%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['40%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['40%']['ivt_2']),
        '50': (tripsect_overlap_dict_out['50%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['50%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['50%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['50%']['ivt_2']), 
        '60': (tripsect_overlap_dict_out['60%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['60%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['60%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['60%']['ivt_2']),
        '70': (tripsect_overlap_dict_out['70%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['70%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['70%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['70%']['ivt_2']), 
        '80': (tripsect_overlap_dict_out['80%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['80%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['80%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['80%']['ivt_2']),
        '95': (tripsect_overlap_dict_out['95%']['ivt_1_2_3'] + 
               tripsect_overlap_dict_out['95%']['ivt_1_2'] + 
               tripsect_overlap_dict_out['95%']['ivt_2_3'] + 
               tripsect_overlap_dict_out['95%']['ivt_2'])}

sim_dict = {}
for i in ['30', '40', '50', '60', '70', '80', '95']:
    sim_dict[int(i)] = {'ivt_1_2':d_sim((tripsect_overlap_dict_out[f"{i}%"]['ivt_1_2_3'] +
                                    tripsect_overlap_dict_out[f"{i}%"]['ivt_1_2']), 
                                   ivt1[i], 
                                   ivt2[i]),
                   'ivt_2_3':d_sim((tripsect_overlap_dict_out[f"{i}%"]['ivt_1_2_3'] + 
                                    tripsect_overlap_dict_out[f"{i}%"]['ivt_2_3']), 
                                   ivt2[i], 
                                   ivt3[i]),
                   'ivt_1_3':d_sim((tripsect_overlap_dict_out[f"{i}%"]['ivt_1_2_3'] + 
                                    tripsect_overlap_dict_out[f"{i}%"]['ivt_1_3']), 
                                   ivt1[i], 
                                   ivt3[i])}
    
df_sim = pd.DataFrame.from_dict(sim_dict)
df_sim = df_sim.T
df_sim.columns = ['IVT 1 & 2 Simmilarity', 'IVT 2 & 3 Simmilarity', 'IVT 1 & 3 Simmilarity']
df_sim.to_csv('Biological_rep_dice_sim.csv')

sns.set(rc={"figure.figsize":(12, 12)}) #width=3, #height=4
sns.set(font_scale=1)
x = sns.lineplot(data=df_sim)
x.set(title='Dice Similarity of SNVs Found in Locations Common to All HeLa Replicates', 
      xlabel='Percent SNV Occurence Threshold', 
      ylabel='Dice Similarity Score', ylim=[0,1.0])

x.get_figure().savefig('NOFILTERDiceSimTripsectionSNV.pdf')