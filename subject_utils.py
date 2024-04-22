import numpy as np
import hashlib

class DNASubject():
    
    def __init__(self, chm, chm_length_morgans, chm_length_snps, maternal, paternal, name=None, sex=None, birthloc=None):
          
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
        
        self.order = sorted(maternal.keys())
        self.maternal = maternal
        self.paternal = paternal  

        # assert all sequences have same length

        if self.name is None:
            maternal_str = "".join(str(x) for x in self.maternal['snps'].tolist())
            paternal_str = "".join(str(x) for x in self.maternal['snps'].tolist())
            all_str = (maternal_str + paternal_str).encode('utf-8')
            self.name = hashlib.md5(all_str).hexdigest()
        
        self.sex = sex if sex is not None else np.random.randint(2, size=1)

        self.birthloc = birthloc if birthloc is not None else (0,0) # default coordinate
        
    def admix(self, breakpoint_probability=None):

        """
        create an admixed haploid from the paternal and maternal sequences.
        Returns:
        haploid_returns: dict with same keys as self.maternal and self.paternal
        """
        num_crossovers = int(np.random.poisson(self.chm_length_morgans))
        
        # debug...
        # num_crossovers=3
        
        #print(num_crossovers)
        haploid_returns = {}
        for key in self.order:
            haploid_returns[key] = np.zeros_like(self.maternal[key])
        
        # edge case of no breaking points.
        if num_crossovers == 0:
            haploid_returns = {}
            select = self.maternal if np.random.rand()>=0.5 else self.paternal
            for key in self.order:
                haploid_returns[key] = select[key].copy()
                
        
        # now, we have to choose crossover points and return 
        # the appropriate admixed sequence.
        
        # the crossover points must be chosen based on the 
        # genetic map and not the snps we have.
        # TODO: gotta make sure to pick the closest point in the snps we have.
        # For now, just uniform over snp length.
        
        else:
            
            breakpoints = np.random.choice(np.arange(1,self.chm_length_snps), 
                                           size=num_crossovers, 
                                           replace=False, 
                                           p=breakpoint_probability)
            breakpoints = np.sort(breakpoints)
            breakpoints = np.concatenate(([0],breakpoints,[self.chm_length_snps]))
            
            #print(breakpoints)
            # select paternal or maternal randomly and apply crossovers.
            choice = np.random.rand()>=0.5
            select = self.maternal if choice else self.paternal
            cross = self.paternal if choice else self.maternal
            for key in self.order:
                haploid_returns[key] = select[key].copy()
                for i in range(1,len(breakpoints)-1,2):
                    begin = breakpoints[i]
                    end = breakpoints[i+1]
                    haploid_returns[key][begin:end] = cross[key][begin:end].copy()
    
        return haploid_returns

    def __repr__(self):
        return self.name

def create_new_subject(subject1, subject2, gen, neo4j_driver, breakpoint_probability=None):

    maternal = subject1.admix(breakpoint_probability)
    paternal = subject2.admix(breakpoint_probability)

    ## Choose birthplace randomly (paternal or maternal) right now.
    choice = np.random.rand()>=0.5
    birthloc = subject1.birthloc if choice else subject2.birthloc
    
    assert subject1.chm == subject2.chm, "Wrong chromosomes being passed!!!"
    
    chm = subject1.chm
    chm_length_morgans = subject1.chm_length_morgans
    chm_length_snps = subject1.chm_length_snps

    new_subject = DNASubject(chm, chm_length_morgans, chm_length_snps, maternal, paternal, birthloc)
    
    node = Node('DNASubject', name=new_subject.name, gen=gen, sex=new_subject.sex)
    neo4j_driver.create(node)

    return new_subject

def distance(p1, p2):
    R = 6371.0 # Radius of the Earth in km
    
    lat1, lon1 = np.radians(p1[0]), np.radians(p1[1])
    lat2, lon2 = np.radians(p2[0]), np.radians(p2[1])
    
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * atan2(np.sqrt(a), np.sqrt(1 - a))
    
    distance = R * c
    return distance

def gaussian_probability(distance, mean, std):
    exponent = -((distance - mean) ** 2) / (2 * std ** 2)
    probability = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(exponent)
    return probability

def mating_probability(subject1, subject2):
    mean_distance = 100  # Mean distance for the Gaussian distribution (adjust as needed)
    std_dev = 1000  # Standard deviation for the Gaussian distribution (adjust as needed)
    distance = distance_between_points(subject1.birthloc, subject2.birthloc)
    probability = gaussian_probability(distance, mean_distance, std_dev)
    return np.asarray(probability)