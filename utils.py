import numpy as np
import pandas as pd
import scipy.interpolate

def sift_vcf_thru_genetic_map(genetic_map, vcf_data):

    """
    get chromosome length in morgans from genetic map.
    Assumes genetic_map is sorted.
    """

    # getting chm number...
    chm = list(set(vcf_data["variants/CHROM"]))
    if len(chm) == 1: chm = chm[0]
    else: raise Exception("Multiple chromosomes in this file!!!")
    chm = chm.lstrip("chr") # in some reference files we have this instead of just chr number. 22 or chr22

    # read in genetic map and subset to chm number.
    genetic_df = pd.read_csv(genetic_map, delimiter="\t", header=None, comment="#", dtype=str)
    genetic_df.columns = ["chm", "pos", "cM"]
    genetic_chm = genetic_df[genetic_df["chm"] == chm]

    if len(genetic_chm) == 0:
        genetic_chm = genetic_df[genetic_df["chm"] == "chr" + chm] # sometimes it is called chr22 instead of 22

    genetic_chm = genetic_chm.astype({"chm":str, "pos":int, "cM":float})

    # get length of chm.
    chm_length_morgans = max(genetic_chm["cM"])/100.0

    # get snp info - snps in the vcf file and their cm values.
    # then compute per position probability of being a breakpoint.
    # requires some interpolation and finding closest positions.
    """
    # 1: Minimum in a sorted array approach and implemented inside admix().
        - O(logn) every call to admix. Note that admix() takes O(n) anyway.
    # 2: Find probabilities using span. - One time computation.

    """
    # This adds 0 overhead to code runtime.
    # get interpolated values of all reference snp positions
    genomic_intervals = scipy.interpolate.interp1d(x=genetic_chm["pos"].to_numpy(), y=genetic_chm["cM"].to_numpy(),fill_value="extrapolate")
    print(genomic_intervals)
    genomic_intervals = genomic_intervals(vcf_data["variants/POS"])
    lengths = genomic_intervals[1:] - genomic_intervals[0:-1]
    bp = lengths / lengths.sum()

    genetic_map_data = {}
    genetic_map_data["chm"] = chm
    genetic_map_data["chm_length_snps"] = vcf_data["calldata/GT"].shape[0]
    genetic_map_data["chm_length_morgans"] = chm_length_morgans
    genetic_map_data["breakpoint_probability"] = bp

    return genetic_map_data

def get_sample_map_data(list_of_samples_pop, vcf_data):
    
    """
    Inputs:
    list_of_samples_pop: tab delimited file with sample, population and no header.
    vcf_data: allel.read_vcf output. It is the reference vcf file information.
    
    Returns:
    sample_map_data: dataframe with sample, population, population_code and index in vcf_data referecnce.
    
    """
    
    # reading sample map
    # sample_map_data = pd.read_csv(sample_map, delimiter="\t",header=None,comment="#")
    # print(sample_map)
    list_of_samples_pop.columns = ["sample", "population"]

    # creating ancestry map into integers from strings
    # id is based on order in sample_map file.
    ancestry_map = {}
    curr = 0
    for i in list_of_samples_pop["population"]:
        if i in ancestry_map.keys():
            continue
        else:
            ancestry_map[i] = curr
            curr += 1
    print("Ancestry map", ancestry_map)
    list_of_samples_pop["population_code"] = np.vectorize(ancestry_map.get)(list_of_samples_pop["population"])

    # getting index of samples in the reference files

    b = vcf_data["samples"]
    a = np.array(list(list_of_samples_pop["sample"]))

    sorter = np.argsort(b)
    indices = sorter[np.searchsorted(b, a, sorter=sorter)]
    list_of_samples_pop["index_in_reference"] = indices
    
    return list_of_samples_pop