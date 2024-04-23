from tqdm import tqdm
import numpy as np
from py2neo import Graph, Node, Relationship, NodeMatcher
from neo4j import GraphDatabase
 
from subject_utils import DNASubject, create_new_subject, mating_probability, softmax, analytical_ancestry


def build_founders(vcf_data, genetic_map, sample_map, neo4j_driver, out_dir, population_mapper):
    
    # Information about snps, chm, lengths from reference, genetic map.
    chm = genetic_map["chm"]
    chm_length_morgans = genetic_map["chm_length_morgans"]
    chm_length_snps = genetic_map["chm_length_snps"]

    # building founders
    founders = []
    index = sample_map.index
    #index = index[:11] #####!!!!
    index = sample_map[sample_map.population =='EUR'].sample(n=5).index.append(sample_map[sample_map.population =='OCE'].sample(n=5).index).append(sample_map[sample_map.population =='EAS'].sample(n=5).index)
    
    print(index)
    
    print(sample_map.loc[index,:])
    for i in tqdm(index):
        
        name = sample_map.loc[i, "sample"]

        # when creating maternal, paternal make sure it has same keys

        maternal = {}
        paternal = {}

        # let us use the first for maternal in the vcf file...
        """
        maternal["snps"] = vcf_data["calldata/GT"][:,i,0]
        paternal["snps"] = vcf_data["calldata/GT"][:,i,1]
        """

        # single ancestry assumption.
        maternal["anc"] = np.array([sample_map.loc[i,'population_code']]*chm_length_snps)
        paternal["anc"] = np.array([sample_map.loc[i,'population_code']]*chm_length_snps)
        
        
        # any more info like coordinates, prs can be added here.
        subject = DNASubject(
            chm=chm, 
            chm_length_morgans=chm_length_morgans, 
            chm_length_snps=chm_length_snps, 
            maternal=maternal, 
            paternal=paternal, 
            birthloc=(sample_map.loc[i, 'latitude'], sample_map.loc[i, 'longitude'])
        )
        founders.append(subject)
        
        ancestry_percentages = analytical_ancestry(maternal["anc"],paternal["anc"],population_mapper)
        
        node = Node(
            'DNASubject_Gen0', 
            name=subject.name, 
            gen=0, 
            sex=subject.sex, 
            progenitor1='root', 
            progenitor2='root',
            **ancestry_percentages
        )
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
    population_mapper,
    panmixia_factor=0.000001,
    seed=42,
):
    print("Building founders")
    # That's generation 0
    founders = build_founders(vcf_data, genetic_map, sample_map, neo4j_driver, out_dir, population_mapper)

    print("Simulating...")
    np.random.seed(seed)
    
    population = {}
    population[0] = founders

    for gen in range(0, num_gens):
        if gen == 0: # nothing to do, it's founders!
            continue

        print("Simulating generation ", gen)
        current_gen_subjects = []

        num_progenitors_gen = len(range(len(population[gen-1])))
        
        print(range(num_samples_per_gen))

        for i in tqdm(range(num_samples_per_gen)):
            # Select parents.

            progenitor1_idx = np.random.choice(num_progenitors_gen)
            progenitor1 = population[gen-1][progenitor1_idx]

            probs_mating = np.zeros(num_progenitors_gen)

            for progenitor2_idx in range(num_progenitors_gen):
                if progenitor1_idx != progenitor2_idx: 

                    # Compute mating probability between progenitor1 and progenitor2.
                    # As of now, it is based on distance. 
                    progenitor2 = population[gen-1][progenitor2_idx]
                    probs_mating[progenitor2_idx] = mating_probability(progenitor1, progenitor2)
            
            probs_mating = np.asarray(probs_mating)
            probs_mating = softmax(probs_mating, temperature=panmixia_factor)
            print('After softmax:', probs_mating)
            # Redistribute probability mass form progenitor1.
            probs_mating += (probs_mating[progenitor1_idx] / (len(probs_mating) - 1))
            probs_mating[progenitor1_idx] = 0

             # Choose mate for progenitor1.
            print('After mass redistribution:', probs_mating)
            print(sum(probs_mating))
            progenitor2_idx = np.random.choice(range(num_progenitors_gen), p=probs_mating) 
            progenitor2 =  population[gen-1][progenitor2_idx]
    
            # non-recursive simulation
            new_subject = create_new_subject(progenitor1, progenitor2, gen, population_mapper, neo4j_driver, breakpoint_probability=None) ###!!!
            current_gen_subjects.append(new_subject)

        population[gen] = current_gen_subjects

    return population

    #if out_root == None:
    #    return dataset # useful when we want to create dataset iterators.
    print("Writing output")
    write_output(out_root, dataset)