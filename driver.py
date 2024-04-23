import os
os.chdir(os.path.join(os.environ.get('HOME'),'genkins'))
os.environ["SRC_PATH"]="/home/users/geleta/genkins"
os.environ["DATA_PATH"]="/scratch/groups/cdbustam/rita/1000G"
os.environ["OUT_PATH"]= os.path.join(os.environ.get('DATA_PATH'), "simulated")

import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

from py2neo import Graph, Node, Relationship, NodeMatcher
from neo4j import GraphDatabase
from simulators import build_trees
from utils import sift_vcf_thru_genetic_map

# Connect to neo4j database.
# https://neo4j.com/docs/api/python-driver/current/api.html#uri-ref
print(os.environ.get('NEO4J_URI'), (os.environ.get('NEO4J_USER'),os.environ.get('NEO4J_PASS')))
graph = Graph(os.environ.get('NEO4J_URI'), auth=(os.environ.get('NEO4J_USER'),os.environ.get('NEO4J_PASS')))

# Load all genetic files needed for simulation.

"""
## Load VCF data.
print('Loading VCF data...')

chr22_vcf_path = os.path.join(os.environ.get('DATA_PATH'), 'ref_final_beagle_phased_1kg_hgdp_sgdp_chr22_hg19.vcf')
chr22_vcf_data = allel.read_vcf(chr22_vcf_path)

## Load genetic map.
print('Loading genetic map...')
genetic_map_path = os.path.join(os.environ.get('DATA_PATH'), 'allchrs.b37.gmap')
chm22_genetic_map = sift_vcf_thru_genetic_map(genetic_map=genetic_map_path, vcf_data=chr22_vcf_data)
"""
chm22_genetic_map = np.load(os.path.join(os.environ.get('DATA_PATH'), 'chm22_genetic_map.npy'), allow_pickle=True).item()

chr22_vcf_data = None

## Load sample map.
print('Loading sample map...')
reference_panel_path = os.path.join(os.environ.get('DATA_PATH'), 'reference_panel_metadata.tsv')
reference_sample_map = pd.read_csv(reference_panel_path, sep="\t")
# Remove unnecessary columns.
reference_sample_map = reference_sample_map.drop(columns=['Population code','Source','Region','Sample Alias','Country', 'Town'])
# Filter to single ancestry.
reference_sample_map = reference_sample_map[reference_sample_map["Single_Ancestry"] == 1]
sample_map = reference_sample_map[['Sample', 'Superpopulation code', 'Latitude', 'Longitude']]
sample_map.columns = ['sample', 'population', 'latitude', 'longitude']
sample_map['population_code'] = sample_map['population'].astype('category')
sample_map['population_code'] = sample_map.population_code.cat.codes
population_mapper = dict(zip(sample_map.population_code, sample_map.population))
print('Population mapper', population_mapper)

#print(sample_map.head())


# Check that we are starting from scratch.
nodes = NodeMatcher(graph)
if len(nodes) > 0: 
    answer = input('Graph DB is not empty!!! Are you sure you want to delete the graph?')
    if answer.lower() in ["y","yes"]:
        query = "MATCH (n) DETACH DELETE n"
        graph.run(query)  
    else:
        exit(1)
    

# Start building trees.
print('Start building trees...')
build_trees(
    vcf_data=chr22_vcf_data, 
    sample_map=sample_map, 
    genetic_map=chm22_genetic_map, 
    num_gens=5,
    num_samples_per_gen=5,
    neo4j_driver=graph,
    out_dir=os.environ.get('OUT_PATH'),
    population_mapper=population_mapper,
    panmixia_factor=0.000001, ## the bigger the panmixia_factor, the higher the panmixia in the population
)