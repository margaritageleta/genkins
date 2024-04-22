scp geleta@galangal.stanford.edu:/home/projects/world_wide_references/reference_panel_metadata.tsv /Users/geleta/Downloads
scp /Users/geleta/Downloads/reference_panel_metadata.tsv geleta@login.sherlock.stanford.edu:/scratch/groups/cdbustam/rita/1000G

scp geleta@galangal.stanford.edu:/home/database/maps/rfmix/allchrs.b37.gmap /Users/geleta/Downloads
scp /Users/geleta/Downloads/allchrs.b37.gmap geleta@login.sherlock.stanford.edu:/scratch/groups/cdbustam/rita/1000G

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	scp geleta@galangal.stanford.edu:/home/projects/world_wide_references/ref_final_$i/ref_final_beagle_phased_1kg_hgdp_sgdp_chr$i\_hg19.vcf.gz /Users/geleta/Downloads
    scp /Users/geleta/Downloads/ref_final_beagle_phased_1kg_hgdp_sgdp_chr$i\_hg19.vcf.gz geleta@login.sherlock.stanford.edu:/scratch/groups/cdbustam/rita/1000G
done


