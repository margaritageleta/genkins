import numpy as np
import hashlib
from py2neo import Graph, Node, Relationship, NodeMatcher
from neo4j import GraphDatabase

class DNASubject():
    
    def __init__(
        self, 
        chm, 
        chm_length_morgans, 
        chm_length_snps, 
        maternal, 
        paternal, 
        name=None, 
        sex=None, 
        birthloc=None,
        marriage_organization='agamous',   ## default
        marriage_composition='monogamous', ## default
        cousin_marriages='any_second_cousins_permitted', ## default
    ):
          
        """
        Inputs:
        chm: chm number.
        chm_length_morgans: chromosome length in morgans found using the genetic map.
        chm_length_snps: chromosome length in avalaible snps based on reference file.
        
        maternal, paternal are dictionaries with the keys: snps, anc, etc... 
        """ 
        
        # chm related information
        self.chm = chm
        self.chm_length_morgans = chm_length_morgans
        self.chm_length_snps = chm_length_snps
        
        assert maternal.keys() == paternal.keys(), "Key mismatch!!!"
        
        self.haploid_keys = sorted(maternal.keys())
        self.maternal = maternal
        self.paternal = paternal  

        # assert all sequences have same length

        #if name is None:
        #    if len(self.maternal) > 0: 
        #        maternal_str = "".join(str(x) for x in self.maternal['snps'].tolist())
        #        paternal_str = "".join(str(x) for x in self.maternal['snps'].tolist())
        #        all_str = (maternal_str + paternal_str).encode('utf-8')
        #    else:
        #        all_str = (str(int(np.random.rand()*1e6))).encode('utf-8')
        #    self.name = hashlib.md5(all_str).hexdigest()
        #else: 
        all_str = str(int(np.random.rand()*1e6))
        all_str = (str(name) + all_str).encode('utf-8')
        self.name = hashlib.md5(all_str).hexdigest() #!!!!!
        
        self.sex = int(sex) if sex is not None else int(np.random.rand()>=0.5)

        self.birthloc = birthloc if birthloc is not None else (0,0) # default coordinate
        self.marriage_organization = marriage_organization if marriage_organization is not None else 'agamous'
        self.marriage_composition = marriage_composition if marriage_composition is not None else 'monogamous'
        self.cousin_marriages = cousin_marriages if cousin_marriages is not None else 'any_second_cousins_permitted'
        
    def admix(self, breakpoint_probability=None):

        """
        create an admixed haploid from the paternal and maternal sequences.
        Returns:
        haploid_returns: dict with same keys as self.maternal and self.paternal
        """
        num_crossovers = int(np.random.poisson(self.chm_length_morgans))
        print('NUM OF CROSSOVERS!', num_crossovers)
        
        # debug...
        # num_crossovers=3
        
        haploid_returns = {}
        
        # edge case of no breaking points.
        if num_crossovers == 0:
            select = self.maternal if np.random.rand()>=0.5 else self.paternal
            for key in self.haploid_keys:
                haploid_returns[key] = select[key].copy()
                
        
        # now, we have to choose crossover points and return 
        # the appropriate admixed sequence.
        
        # the crossover points must be chosen based on the 
        # genetic map and not the snps we have.
        # TODO: gotta make sure to pick the closest point in the snps we have.
        # For now, just uniform over snp length.
        
        else:
            breakpoints = np.random.choice(
                np.arange(1,self.chm_length_snps), 
                size=num_crossovers, 
                replace=False, 
                p=breakpoint_probability
            )
            breakpoints = np.sort(breakpoints)
            breakpoints = np.concatenate(([0],breakpoints,[self.chm_length_snps]))
            
            #print(breakpoints)
            # select paternal or maternal randomly and apply crossovers.
            choice = np.random.rand() >= 0.5
            select = self.maternal if choice else self.paternal
            cross = self.paternal if choice else self.maternal
            for key in self.haploid_keys:
                haploid_returns[key] = select[key].copy()
                for i in range(1, len(breakpoints) - 1, 2):
                    begin = breakpoints[i]
                    end = breakpoints[i + 1]
                    haploid_returns[key][begin:end] = cross[key][begin:end].copy()
    
        return haploid_returns

    def __repr__(self):
        return f"""
        Subject {self.name} (chm{self.chm})
        ------------------
        Sex: {self.sex}
        Birth location: {self.birthloc}
        Marriage organization: {self.marriage_organization}
        Marriage composition: {self.marriage_composition}
        Cousin marriages: {self.cousin_marriages}
        """

def create_new_subject(subject1, subject2, gen, population_mapper, neo4j_driver, breakpoint_probability=None):

    maternal = subject1.admix(breakpoint_probability)
    paternal = subject2.admix(breakpoint_probability)

    ## Choose birthplace randomly (paternal or maternal) right now.
    choice = np.random.rand()>=0.5
    birthloc = subject1.birthloc if choice else subject2.birthloc
    sex = bool(np.random.rand()>=0.5)
    
    assert subject1.chm == subject2.chm, "Wrong chromosomes being passed!!!"
    
    chm = subject1.chm
    chm_length_morgans = subject1.chm_length_morgans
    chm_length_snps = subject1.chm_length_snps

    new_subject = DNASubject(
        chm=chm, 
        chm_length_morgans=chm_length_morgans, 
        chm_length_snps=chm_length_snps, 
        maternal=maternal, 
        paternal=paternal, 
        birthloc=birthloc,
        name=gen,
        marriage_organization=subject1.marriage_organization if sex else subject2.marriage_organization,
        marriage_composition=subject1.marriage_composition if sex else subject2.marriage_composition,
        cousin_marriages=subject1.cousin_marriages if sex else subject2.cousin_marriages,
    )
    
    # Add child to graph.
    ancestry_percentages = analytical_ancestry(maternal["anc"],paternal["anc"],population_mapper)
    
    node = Node(
        f'DNASubject', 
        name=new_subject.name, 
        gen=gen, 
        sex=new_subject.sex, 
        progenitor1=subject1.name, 
        progenitor2=subject2.name, 
        **ancestry_percentages,
        marriage_organization=new_subject.marriage_organization,
        marriage_composition=new_subject.marriage_composition,
        cousin_marriages=new_subject.cousin_marriages,
    )
    neo4j_driver.create(node)
    
    # Connect child to progenitors.
    nodes = NodeMatcher(neo4j_driver)
    progenitor1 = nodes.match(f'DNASubject', name=subject1.name, gen=gen-1).first()
    progenitor2 = nodes.match(f'DNASubject', name=subject2.name, gen=gen-1).first()
    neo4j_driver.create(Relationship(progenitor1, 'PARENT', node))
    neo4j_driver.create(Relationship(progenitor2, 'PARENT', node))

    return new_subject

def softmax(x, temperature=0.000001):
    """Compute softmax values for each element of input array with temperature."""
    max_x = np.max(x)
    exp_x = np.exp((x - max_x) / temperature)
    sum_exp_x = np.sum(exp_x)
    return exp_x / sum_exp_x

def distance(p1, p2):
    R = 6371.0 # Radius of the Earth in km
    
    lat1, lon1 = np.radians(p1[0]), np.radians(p1[1])
    lat2, lon2 = np.radians(p2[0]), np.radians(p2[1])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    return R * c

def gaussian_probability(distance, mean, std):
    exponent = -((distance - mean) ** 2) / (2 * std ** 2)
    probability = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(exponent)
    return probability

def mating_probability(subject1, subject2):
    mean_distance = 100  # Mean distance for the Gaussian distribution (adjust as needed)
    std_dev = 1000  # Standard deviation for the Gaussian distribution (adjust as needed)
    d = distance(subject1.birthloc, subject2.birthloc)
    probability = gaussian_probability(d, mean_distance, std_dev)
    return probability

def analytical_ancestry(maternal, paternal, population_mapper):
    full_strand = np.concatenate([maternal, paternal])
    full_strand_len = len(full_strand)
    ancestry = {}
    for i, anc in population_mapper.items():
        ancestry[f'ancestry_{anc}'] = len(np.where(full_strand == i)[0]) / full_strand_len * 100
    return ancestry