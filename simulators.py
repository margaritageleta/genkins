from tqdm import tqdm
import numpy as np
from py2neo import Graph, Node, Relationship, NodeMatcher
from neo4j import GraphDatabase
 
from subject_utils import DNASubject, create_new_subject, mating_probability

def build_founders(vcf_data, genetic_map, sample_map, neo4j_driver):
    
    # Information about snps, chm, lengths from reference, genetic map.
    chm = genetic_map["chm"]
    chm_length_morgans = genetic_map["chm_length_morgans"]
    chm_length_snps = genetic_map["chm_length_snps"]

    # building founders
    founders = []

    for i in tqdm(range(sample_map.shape[0])):

        index = sample_map.loc[i, "index_in_reference"]
        name = sample_map.loc[i, "sample"]

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
        subject = DNASubject(chm, chm_length_morgans, chm_length_snps, maternal, paternal)
        founders.append(subject)
        
        node = Node('DNASubject', name=subject.name, gen=0, sex=subject.sex)
        neo4j_driver.create(node)

    

    return founders

def build_trees(
    vcf_data, 
    sample_map, 
    genetic_map, 
    num_gens,
    num_samples_per_gen, 
    neo4j_driver,
    out_dir, 
    seed=42
):
    print("Building founders")
    # That's generation 0
    founders = build_founders(vcf_data, genetic_map, sample_map, neo4j_driver, out_dir)

    print("Simulating...")
    np.random.seed(seed)
    
    population = {}
    population[0] = founders

    for gen in range(0, num_gens):
        if gen == 0: # nothing to do, it's founders!
            continue

        print("Simulating generation ", gen)
        current_gen_subjects = []

        num_progenitors_gen = len(range(population[gen-1]))

        for i in range(num_samples_per_gen):
            # Select parents.

            progenitor1_idx = random.choice(num_progenitors_gen)
            progenitor1 = population[gen-1][progenitor1_idx]

            probs_mating = np.zeros(num_progenitors_gen)

            for progenitor2_idx in range(num_progenitors_gen):
                if prorgenitor1_idx != progenitor2_idx: 

                    # Compute mating probability between progenitor1 and progenitor2.
                    # As of now, it is based on distance. 
                    progenitor2 = population[gen-1][progenitor2_idx]
                    probs_mating[progenitor2_idx] = mating_probability(progenitor1, progenitor2)
            
            # Redistribute probability mass form progenitor1.
            probs_mating += (probs_mating[progenitor1_idx] / (len(probs_mating) - 1))
            probs_mating[progenitor1_idx] = 0

             # Choose mate for progenitor1.
            print(probs_mating)
            progenitor2_idx = np.random.choice(range(len(num_progenitors_gen)), p=probs_mating) 
            progenitor2 =  population[gen-1][progenitor2_idx]
    
            # non-recursive simulation
            new_subject = create_subject(progenitor1, progenitor2, gen, breakpoint_probability)
            current_gen_subjects.append(new_subject)

        population[gen] = current_gen_subjects

    return population

    #if out_root == None:
    #    return dataset # useful when we want to create dataset iterators.
    print("Writing output")
    write_output(out_root, dataset)