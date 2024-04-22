from tqdm import tqdm
import numpy as np
 
from subject_utils import DNASubject, create_new_subject

def build_founders(vcf_data, genetic_map, sample_map):
    
    # information about snps, chm, lengths from reference, genetic map.
    chm = genetic_map_data["chm"]
    chm_length_morgans = genetic_map_data["chm_length_morgans"]
    chm_length_snps = genetic_map_data["chm_length_snps"]

    # building founders
    founders = []

    for i in tqdm(range(sample_map.shape[0])):

        # first get the index of this sample in the vcf_data.
        # if not there, skip and print to log.

        index = sample_map.loc[i,"index_in_reference"]
        name = sample_map.loc[i,"sample"]

        # when creating maternal, paternal make sure it has same keys

        maternal = {}
        paternal = {}

        # let us use the first for maternal in the vcf file...
        maternal["snps"] = vcf_data["calldata/GT"][:,index,0]
        paternal["snps"] = vcf_data["calldata/GT"][:,index,1]

        # single ancestry assumption.
        maternal["anc"] = np.array([sample_map.loc[i,"population_code"]]*chm_length_snps)
        paternal["anc"] = np.array([sample_map.loc[i,"population_code"]]*chm_length_snps)


        # any more info like coordinates, prs can be added here.

        subject = DNASubject(chm, chm_length_morgans, chm_length_snps, maternal, paternal, name)
        founders.append(subject)

    return founders

def simulate(vcf_data, sample_map, genetic_map, num_gens,
             num_samples_per_gen, out_root=None, seed=42):

    print("Building founders")
    founders = build_founders(vcf_data, genetic_map, sample_map)
    return founders

    print("Simulating...")
    dataset = build_tree(
        founders=founders,
        num_gens=num_gens,
        num_samples_per_gen=num_samples_per_gen,
        breakpoint_probability=genetic_map["breakpoint_probability"],
        seed=seed,
    )
    if out_root == None:
        return dataset # useful when we want to create dataset iterators.
    print("Writing output")
    write_output(out_root, dataset)