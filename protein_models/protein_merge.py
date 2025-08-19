#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:
prefix = r'protein_dataset/'
reproduced_dict = {'esm1b_t33_650M_UR50S' : {'name' : 'ESM1b',
                                             'score_path' : prefix + 'ESM1b_score.csv'
                                            },
                   'esm1v_t33_650M_UR90S_1' : {'name' : 'ESM1v-1',
                                               'score_path' : prefix + 'ESM1v_score.csv'
                                              },
                   'esm2_t33_650M_UR50D' : {'name' : 'ESM2',
                                            'score_path' : prefix + 'ESM2_score.csv'
                                           }
                  }

output_path = r'./protein_dataset/ESM_score.csv'


# In[3]:



# In[5]:


merge_df = None
for model_pretrained, model_dict in reproduced_dict.items():
    score_df = pd.read_csv(model_dict['score_path'])

    if 'esm_score' in score_df.columns:
        score_df = score_df.rename(columns={'esm_score': model_dict['name'] + '_score'})

    if merge_df is None:
        merge_df = score_df
    else:
        merge_df = pd.merge(merge_df, score_df, on=['seq_id', 'mut_name'], how='inner')

# merge_df = reduce(lambda left, right: pd.merge(left, right, on=['seq_id', 'mut_name'], how='inner'), model_df_list)
merge_df.to_csv(output_path, index=False)


# In[ ]:




