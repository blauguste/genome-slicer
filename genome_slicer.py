from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
from collections import defaultdict

def find_genbank_feature(record, locus_tag):
    for feature in record.features:
        if feature.type == 'gene':
            if feature.qualifiers['locus_tag'][0] == locus_tag:
                corresponding_gb_feature = feature
    return corresponding_gb_feature

def make_IGR_feature(start_locus, end_locus, igr_feature_type):
    igr_feature_location = FeatureLocation(start_locus, end_locus, strand=None)
    igr_feature = SeqFeature(igr_feature_location, type=igr_feature_type, strand=None)
    return igr_feature

def get_ss_IGRs(gene_list, min_igr_length, igr_feature_type):
    ss_IGR_set = set()
    # Sort list of genes by start position
    gene_list.sort(key=lambda tup: tup[1])
    base_gene_list_position = 0
    comparison_gene_list_position = 1
    # Compare the right boundary of the "first" gene with the left boundary of the subsequent gene
    for i, gene_record in enumerate(gene_list[comparison_gene_list_position:]):
        base_gene_right = gene_list[base_gene_list_position][2]
        comparison_gene_left = gene_list[comparison_gene_list_position][1]
        comparison_gene_right = gene_list[comparison_gene_list_position][2]
        base_gene_name = gene_list[base_gene_list_position][0]
        comparison_gene_name = gene_list[comparison_gene_list_position][0]
        # If there's a gap of desired size between the two genes...
        if comparison_gene_left - base_gene_right >= min_igr_length:
            # Create new features
            igr_feature = make_IGR_feature(base_gene_right, comparison_gene_left, igr_feature_type)
            igr_feature.qualifiers["locus_tag"] = base_gene_name + "_IGR_" + comparison_gene_name
            ss_IGR_set.add(igr_feature)
            # Reset the base gene and comparison position (move to next gene)
            base_gene_list_position = comparison_gene_list_position
            comparison_gene_list_position = base_gene_list_position + 1
        # If the right boundary of the comparison gene is greater than that
        # of the base gene, make the comparison gene the new base gene
        elif comparison_gene_right > base_gene_right:
            base_gene_list_position = comparison_gene_list_position
            comparison_gene_list_position = base_gene_list_position + 1
        # If none of these conditions are met, repeat with the next comparison gene
        else:
            comparison_gene_list_position += 1
    return ss_IGR_set

with open('229.opr', 'r') as operon_input:
    operon_reader = csv.reader(operon_input, delimiter = '\t')
    next(operon_reader, None) # ignore the header row
    input_list = []
    for row in operon_reader:
        input_list.append((row[0], row[2]))

record = SeqIO.read('NC_004350_2.gb', 'genbank')
print(record.description, record.id)

# Make a complete list of all the gene features in the genome, store as a set.
complete_feature_set = set()
for feature in record.features:
    if feature.type == 'gene':
        complete_feature_set.add(feature)
print('gene feature list starting length: ', len(complete_feature_set))

# Create a dictionary where the key is operon ID and the value is 
#### a list of gene features contained in that operon
genes_in_operons = set()
operon_dict = defaultdict(list)
for operon_id, locus_tag in input_list:
    operon_feature = find_genbank_feature(record, locus_tag)
    operon_dict[operon_id].append(operon_feature)
    # Keep track of all the gene features that have been added to operons and add them to a set
    genes_in_operons.add(operon_feature)

# Loop through the operon dictionary and create a feature for each operon.
operon_list = []
for operon_id, gene_feature_list in operon_dict.items():
    start = gene_feature_list[0].location.start
    end = gene_feature_list[-1].location.end
    strand = gene_feature_list[0].strand
    new_operon_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type='operon')
    operon_list.append(new_operon_feature)

# Identify all the gene features that were in the original genbank but are NOT among those added to operons
gene_features_not_in_operons = complete_feature_set.difference(genes_in_operons)
print('gene features NOT in operons: ', len(gene_features_not_in_operons))

################################### DEFINE REGION 1: IGRs within operons ###################################

# Generate the gene list needed to feed the ss_IGR finder function (one list per operon)
igrs_within_operons = set()
for operon_id, gene_feature_list in operon_dict.items():
    operon_gene_list = []
    for entry in gene_feature_list:
        start = entry.location.start
        end = entry.location.end
        gene_name = entry.qualifiers["locus_tag"][0]
        strand = entry.strand
        operon_gene_list.append((gene_name, start, end, strand))
    # Search for IGRs within each operon.
    operon_igrs = get_ss_IGRs(operon_gene_list, 1, 'intra_operon_IGR')
    igrs_within_operons.update(operon_igrs)
print('intra-operon IGRs identified: ', len(igrs_within_operons))

############################# DEFINE REGION 2: IGRs between transciptional units ############################

# Create a gene set consisting of 