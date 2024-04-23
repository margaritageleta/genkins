from tqdm import tqdm
import numpy as np
from py2neo import Graph, Node, Relationship, NodeMatcher
from neo4j import GraphDatabase
 
from subject_utils import DNASubject, create_new_subject, mating_probability, softmax, analytical_ancestry


# Check if subject has mates.
def check_current_mates(subject, gen, neo4j_driver):
    nodes = NodeMatcher(neo4j_driver)
    
    # Check when subject is progenitor1.
    offspring1 = nodes.match(f'DNASubject_Gen{gen + 1}', progenitor1=subject.name)
    offspring2 = nodes.match(f'DNASubject_Gen{gen + 1}', progenitor2=subject.name)
    
    mates = {}
    for o in offspring1:
        if o['progenitor2'] in mates:
            mates[o['progenitor2']] += 1
        else: mates[o['progenitor2']] = 1
    for o in offspring2:
        if o['progenitor1'] in mates:
            mates[o['progenitor1']] += 1
        else: mates[o['progenitor1']] = 1
    print(f'Subject {subject.name} has mates: {mates}')
    
    return mates

# Find DNASubject class by name from neo4j:
def get_subject_by_name(name, population):
    num_progenitors_gen = len(range(len(population)))
    subject = None
    for idx in range(num_progenitors_gen):
        if population[idx].name == name:
            subject = population[idx]
    return subject

# If subject is in search of mate:
def find_mate(subject_idx, gen, population, neo4j_driver, panmixia_factor=0.0001):
    
    subject = population[subject_idx]
    num_progenitors_gen = len(range(len(population)))
    
    probs_mating = np.zeros(num_progenitors_gen)

    for mate_idx in range(num_progenitors_gen):
        if subject_idx != mate_idx: 

            # Compute mating probability between progenitor1 and progenitor2.
            # As of now, it is based on distance. 
            mate = population[mate_idx]
            probs_mating[mate_idx] = mating_probability(subject, mate)

    probs_mating = np.asarray(probs_mating)
    probs_mating = softmax(probs_mating, temperature=panmixia_factor)
    print('After softmax:', probs_mating)
    
    # Redistribute probability mass from subject.
    probs_mating += (probs_mating[subject_idx] / (len(probs_mating) - 1))
    probs_mating[subject_idx] = 0

     # Choose mate for subject.
    print('After mass redistribution:', probs_mating)
    print(sum(probs_mating))
    mate_idx = np.random.choice(range(num_progenitors_gen), p=probs_mating) 
    mate =  population[mate_idx]
    
    return mate
    


def build_founders(vcf_data, genetic_map, sample_map, neo4j_driver, out_dir, population_mapper):
    
    # Information about snps, chm, lengths from reference, genetic map.
    chm = genetic_map["chm"]
    chm_length_morgans = genetic_map["chm_length_morgans"]
    chm_length_snps = genetic_map["chm_length_snps"]

    # building founders
    founders = []
    sample_map = sample_map[sample_map.marriage_composition == 'limited_polygyny']
    index = sample_map.index
    #index = index[:11] #####!!!!
    index = sample_map[sample_map.population =='AFR'].sample(n=5).index.append(sample_map[sample_map.population =='EAS'].sample(n=5).index)
    
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
            birthloc=(sample_map.loc[i, 'latitude'], sample_map.loc[i, 'longitude']),
            marriage_organization=sample_map.loc[i,'marriage_organization'],
            marriage_composition=sample_map.loc[i,'marriage_composition'],
            cousin_marriages=sample_map.loc[i,'cousin_marriages'],
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
            **ancestry_percentages,
            marriage_organization=subject.marriage_organization,
            marriage_composition=subject.marriage_composition,
            cousin_marriages=subject.cousin_marriages,
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
    nodes = NodeMatcher(neo4j_driver)
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

            ## Select progenitor1.
            progenitor1_idx = np.random.choice(num_progenitors_gen)
            progenitor1 = population[gen-1][progenitor1_idx]
            
            # Check if progenitor1 has mates.
            mates = check_current_mates(progenitor1, gen-1, neo4j_driver) ## generation of progenitor1
            """
            mates = { name_of_mate : #_of_kids_together, ... }
            """
############## CASE 1: If monogamous. #####################################################
            if progenitor1.marriage_composition == 'monogamous':
                # Already has mate:    
                if len(mates.keys()) == 1:
                    progenitor2_name = list(mates.keys())[0]
                    progenitor2 = get_subject_by_name(progenitor2_name, population[gen-1])
                # Has no mate:
                elif len(mates.keys()) == 0: 
                    progenitor2 = find_mate(progenitor1_idx, gen, population[gen-1], neo4j_driver)
                # Monogamous but has lovers??
                else: 
                    raise Exception('I think someone had an affair...')
############## CASE 2: If polygamous. #####################################################               
            elif (progenitor1.marriage_composition == 'limited_polygyny') or (progenitor1.marriage_composition == 'polygyny'):
                if progenitor1.marriage_composition == 'limited_polygyny':
                    polygyny_threshold =  0.85 # probability of starting a new affair is 0.15
                elif progenitor1.marriage_composition == 'polygyny':
                    polygyny_threshold =  0.5  # probability of starting a new affair is 0.5
                # Already has mates: 
                if len(mates.keys()) > 0:
                    start_new_affair = bool(np.random.rand()>=polygyny_threshold) 
                    if start_new_affair:
                        progenitor2 = find_mate(progenitor1_idx, gen, population[gen-1], neo4j_driver)
                    else: # mate with current mates.
                        # weight probability of mate based on # of kids together.
                        mate_options, mate_weigths = mates.keys(), mates.values()
                        mate_options = list(mate_options)
                        mate_weigths = np.asarray(list(mate_weigths))
                        mate_weigths = softmax(mate_weigths, temperature=1)
                        
                        progenitor2_idx = np.random.choice(range(len(mate_options)), p=mate_weigths) 
                        progenitor2_name = mate_options[progenitor2_idx]
                        progenitor2 = get_subject_by_name(progenitor2_name, population[gen-1])
                # Has no mate:
                progenitor2 = find_mate(progenitor1_idx, gen, population[gen-1], neo4j_driver)

############## CASES COVERED: MAKE KIDS ###################################################                         
            new_subject = create_new_subject(progenitor1, progenitor2, gen, population_mapper, neo4j_driver, breakpoint_probability=None) ###!!!
            current_gen_subjects.append(new_subject)

        population[gen] = current_gen_subjects

    return population

    #if out_root == None:
    #    return dataset # useful when we want to create dataset iterators.
    print("Writing output")
    write_output(out_root, dataset)