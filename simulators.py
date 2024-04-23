from tqdm import tqdm
import numpy as np
from py2neo import Graph, Node, Relationship, NodeMatcher
from neo4j import GraphDatabase
 
from subject_utils import DNASubject, create_new_subject, mating_probability, softmax, analytical_ancestry


# Check if subject has mates.
def check_current_mates(subject, gen, neo4j_driver):
    nodes = NodeMatcher(neo4j_driver)
    
    # Check when subject is progenitor1.
    offspring1 = nodes.match(f'DNASubject', progenitor1=subject.name)
    offspring2 = nodes.match(f'DNASubject', progenitor2=subject.name)
    
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

# Check if subjects can mate based on cousin marriage restrictions:
# mate_candidates is a set containing all current candidates.
def prune_cousins(subject, gen, mate_candidates, neo4j_driver):
    nodes = NodeMatcher(neo4j_driver)

    siblings_query = '''
    MATCH (parent1)-[:PARENT]->(subject{{name:"{subject_name}", gen:{generation}}}), (parent2)-[:PARENT]->(sibling{{gen:{generation} }})
    WHERE subject.name <> sibling.name
    AND (subject.progenitor1 = sibling.progenitor1
    OR subject.progenitor1 = sibling.progenitor2
    OR subject.progenitor2 = sibling.progenitor1
    OR subject.progenitor2 = sibling.progenitor2)
    RETURN DISTINCT sibling
    '''
    first_cousins_query = f'''
    
    '''
    
    siblings = [x for x in neo4j_driver.run(siblings_query.format(subject_name=subject.name, generation=gen))]
    siblings = set([sibling['sibling']['name'] for sibling in siblings])
    print(siblings)

    if subject.cousin_marriages == 'any_first_cousins_permitted':
        mate_candidates = mate_candidates.difference(siblings)
        return mate_candidates

    elif subject.cousin_marriages == 'any_second_cousins_permitted':
        mate_candidates = mate_candidates.difference(siblings)
        mate_candidates = mate_candidates.difference(first_cousins)
        return mate_candidates
    
    elif subject.cousin_marriages == 'any_third_cousins_permitted':
        mate_candidates = mate_candidates.difference(siblings)
        mate_candidates = mate_candidates.difference(first_cousins)
        mate_candidates = mate_candidates.difference(second_cousins)
        return mate_candidates
        
        return True # subject and mate are +3rd cousins.

    elif subject.cousin_marriages == 'any_cross_cousins_permitted':
        raise Exception('Not implemented yet.')

    elif subject.cousin_marriages == 'matrilateral_cross_cousins_permitted':
        raise Exception('Not implemented yet.')

    else: # no restrictions...
        # Central thesis of Sigmund Freud and Claude Levi-Strauss: 
        # that exogamous marriage and the establishment of incest taboos
        # are the first acts of the morality and culture which define our species.
        return mate_candidates


# If subject is in search of mate:
def find_mate(subject_idx, gen, population, neo4j_driver, panmixia_factor=0.0001):
    nodes = NodeMatcher(neo4j_driver)
    
    subject = population[subject_idx]
    num_progenitors_gen = len(range(len(population)))
    population_nodes = nodes.match(f'DNASubject', gen=gen).all()
    assert num_progenitors_gen == len(population_nodes), "Mismatch population and graph!!!"
    
    # probs_mating = np.zeros(num_progenitors_gen)

    # In this list we will add all the mate candidates based on the
    # restrictions on cousin marriages and each candidate will
    # have a probability of mating associated based on marriage
    # organization and distance.
    
    mate_candidates = set()
    for mate in population_nodes:
        if subject.name != mate['name']: 
            if subject.sex != mate['sex']:
                mate_candidates.add(mate['name'])
    
    mate_candidates = prune_cousins(subject, gen, mate_candidates, neo4j_driver)
    
    maux = []
    for s in population:
        if s.name in mate_candidates:
            maux.append(s)
    mate_candidates = maux
    if len(mate_candidates) == 0: return None
          
    # Compute mating probability between progenitor1 and progenitor2.
    # As of now, it is based on distance. 
    probs_mating = np.zeros(len(mate_candidates))
    for i, mate in enumerate(mate_candidates):
        probs_mating[i] = mating_probability(subject, mate)

    probs_mating = np.asarray(probs_mating)
    print('Pre softmax:', probs_mating)
    probs_mating = softmax(probs_mating, temperature=panmixia_factor)
    print('After softmax:', probs_mating)

     # Choose mate for subject.
    print('After mass redistribution:', probs_mating)
    print(sum(probs_mating))
    mate_idx = np.random.choice(range(len(mate_candidates)), p=probs_mating) 
    mate =  mate_candidates[mate_idx]
    
    return mate
    


def build_founders(vcf_data, genetic_map, sample_map, neo4j_driver, out_dir, population_mapper):
    
    # Information about snps, chm, lengths from reference, genetic map.
    chm = genetic_map["chm"]
    chm_length_morgans = genetic_map["chm_length_morgans"]
    chm_length_snps = genetic_map["chm_length_snps"]

    # building founders
    founders = []
    sample_map = sample_map[sample_map.cousin_marriages == 'any_first_cousins_permitted']
    sample_map = sample_map[sample_map.marriage_composition == 'monogamous']
    index = sample_map.index
    #index = index[:11] #####!!!!
    #index = sample_map[sample_map.population =='AFR'].sample(n=5).index.append(sample_map[sample_map.population =='EAS'].sample(n=5).index)
    index = sample_map.sample(n=10).index
    
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
            'DNASubject', 
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
                    progenitor2 = find_mate(progenitor1_idx, gen-1, population[gen-1], neo4j_driver)
                # Monogamous but has lovers??
                else: 
                    #raise Exception('I think someone had an affair...')
                    print('I think someone had an affair...')
                    continue
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
                        progenitor2 = find_mate(progenitor1_idx, gen-1, population[gen-1], neo4j_driver)
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
                else: progenitor2 = find_mate(progenitor1_idx, gen-1, population[gen-1], neo4j_driver)

############## CASES COVERED: MAKE KIDS ###################################################  
            if progenitor2 is None:
                print(f'Subject could not find mat. Subject could not find love in life...')
                continue
            new_subject = create_new_subject(progenitor1, progenitor2, gen, population_mapper, neo4j_driver, breakpoint_probability=None) ###!!!
            current_gen_subjects.append(new_subject)

        population[gen] = current_gen_subjects

    return population

    #if out_root == None:
    #    return dataset # useful when we want to create dataset iterators.
    print("Writing output")
    write_output(out_root, dataset)