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


genetic_map_path = os.path.join(os.environ.get('DATA_PATH'), 'allchrs.b37.gmap')
reference_panel_path = os.path.join(os.environ.get('DATA_PATH'), 'reference_panel_metadata.tsv')
chr22_vcf_path = os.path.join(os.environ.get('DATA_PATH'), 'ref_final_beagle_phased_1kg_hgdp_sgdp_chr22_hg19.vcf.gz')

chr22_vcf_data = allel.read_vcf(chr22_vcf_path)
chm22_genetic_map = sift_vcf_thru_genetic_map(genetic_map=genetic_map_path, vcf_data=chr22_vcf_data)

founders = simulate(
    vcf_data=chr22_vcf_data, 
    sample_map=sample_map, 
    genetic_map=chm22_genetic_map, 
    out_root=None,#os.environ.get('OUT_PATH'),
    num_gens=max_gen,
    num_samples_per_gen=num_samples_per_gen, 
    seed=2024
)