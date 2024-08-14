import csv
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

import sys
sys.path.append('/Users/geleta/Developer/genkins')

from dmutils import relabel_country, cluster_islands, cluster_country

interest_vars = [
 'profileid',
 'gender',
 'birth_year',
 'birth_location_city',
 'birth_location_state',
 'birth_location_country',
 'death_year',
 'death_location_city',
 'death_location_state',
 'death_location_country']

print('Reading profiles...')
df = pd.read_csv('../familinx/complete-v2-profiles-anon.txt', sep='\t')
df = df[interest_vars]
df = df.dropna(subset=['profileid', 'gender', 'birth_year', 'death_year', 'birth_location_country'])
print('Processing profiles...')
df['gender'] = pd.Categorical(df['gender'])
df['birth_year'] = pd.to_numeric(df['birth_year'])
df['death_year'] = pd.to_numeric(df['death_year'])
df['age'] = df['death_year'] - df['birth_year']
df = df[df['age'] >= 0]
df = df[df['age'] < 160]
df = df[df['birth_year'] > 1500]
df = df[df['birth_year'] < 2024]

df['birth_location_country'] = df['birth_location_country'].apply(lambda x: relabel_country(x))
df['birth_location_country'] = df['birth_location_country'].apply(lambda x: relabel_country(x))
df['birth_location_country'] = df['birth_location_country'].apply(lambda x: cluster_islands(x))
df['death_location_country'] = df['death_location_country'].apply(lambda x: relabel_country(x))
df['death_location_country'] = df['death_location_country'].apply(lambda x: cluster_islands(x))
df['birth_location_continent'] = df['birth_location_country'].apply(lambda x: cluster_country(x))

gender_counts = df.groupby('birth_location_continent')['gender'].value_counts().unstack(fill_value=0)
gender_counts.columns = ['male_count', 'female_count']
filtered_countries = gender_counts[(gender_counts['male_count'] >= 100) & (gender_counts['female_count'] >= 100)]
df100 = df[df['birth_location_continent'].isin(filtered_countries.index)]

print('Extracting relations...')
profile_ids = df100['profileid'].unique()
profile_ids_set = set(profile_ids)
df100['profileid'] = pd.to_numeric(df100['profileid'])
df100 = df100.set_index('profileid')
edges = dict(zip(profile_ids, [[[],0] for i in range(len(profile_ids))]))
parents_info = {}

with open('../familinx/relations-anon.txt') as csvfile:
    relations_read = 0
    prev_parent_read = -1
    for i, row in enumerate(pbar := tqdm(csv.reader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE))):
        if i == 0: continue
        
        parent, child = row
        parent = int(parent)
        child = int(child)

        if (prev_parent_read == parent) or (parent in profile_ids_set):
            prev_parent_read = parent
            if child in profile_ids_set:
                edges[parent][0].append(child)
                if child in parents_info:
                    if df100.loc[int(parent)].gender == 'female':
                        parents_info[child]['mother'].append(parent)
                    elif df100.loc[int(parent)].gender == 'male':
                        parents_info[child]['father'].append(parent)
                    else: raise Exception('Missing biological sex label.')
                else:
                    parents_info[child] = {
                        'mother': [parent] if df100.loc[int(parent)].gender == 'female' else [],
                        'father': [parent] if df100.loc[int(parent)].gender == 'male' else []
                    }
            edges[parent][1] += 1
            relations_read += 1
        pbar.set_description(f"Processing {i+1} relations ({np.round(min(i+1, 86127511)/86127511 * 100, 2)}%) | Total relations so far: {relations_read} ({np.round(relations_read/(i+1),2 * 100)}%)")

print('Saving relations.')
with open('../familinx/complete-v2-relations-anon.pickle', 'wb') as handle:
    pickle.dump(edges, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('../familinx/complete-v2-parents-info.pickle', 'wb') as handle:
    pickle.dump(parents_info, handle, protocol=pickle.HIGHEST_PROTOCOL)

