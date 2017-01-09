from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
from collections import defaultdict

def feature_to_gene_list(iterable):
    gene_list = []
    for entry in iterable:
        start = entry.location.start
        end = entry.location.end
        gene_name = entry.qualifiers['locus_tag']
        strand = entry.strand
        gene_list.append((gene_name, start, end, strand))
    return gene_list

def find_genbank_feature(record, locus_tag):
    for feature in record.features:
        if feature.type == 'gene':
            if feature.qualifiers['locus_tag'][0] == locus_tag:
                corresponding_gb_feature = feature
    return corresponding_gb_feature

def make_IGR_feature(start_locus, end_locus, igr_feature_type, strand):
    igr_feature_location = FeatureLocation(start_locus, end_locus, strand=strand)
    igr_feature = SeqFeature(igr_feature_location, type=igr_feature_type, strand=strand)
    return igr_feature

def get_ss_IGRs(ref, gene_list, min_igr_length, igr_feature_type, strand, leading_trailing=False):
    ss_IGR_set = set()   
    # Sort list of genes by start position
    gene_list.sort(key=lambda tup: tup[1])
    base_gene_list_position = 0
    comparison_gene_list_position = 1
    if leading_trailing == True:
        # Make the IGR between the origin and the 'first' gene.
        if gene_list[0][1] > 0:
            first_gene_left = gene_list[0][1]
            first_gene_name = gene_list[0][0]
            start_IGR = make_IGR_feature(0, first_gene_left, igr_feature_type, strand)
            start_IGR.qualifiers["locus_tag"] = "start_IGR_" + str(first_gene_name)
            ss_IGR_set.add(start_IGR)
        # Sort the gene list by end position
        gene_list.sort(key=lambda tup: tup[2])
        # Make the the IGR between 'last' gene and 'end' of genome
        if gene_list[-1][1] < len(ref):
            last_gene_right = gene_list[-1][2]
            last_gene_name = gene_list[-1][0]
            end_IGR = make_IGR_feature(last_gene_right, len(ref), igr_feature_type, strand)
            end_IGR.qualifiers["locus_tag"] = str(last_gene_name) + "_IGR_end"
            ss_IGR_set.add(end_IGR)
        # Sort list of genes by start position
        gene_list.sort(key=lambda tup: tup[1])
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
            igr_feature = make_IGR_feature(base_gene_right, comparison_gene_left, igr_feature_type, strand)
            igr_feature.qualifiers["locus_tag"] = str(base_gene_name) + "_IGR_" + str(comparison_gene_name)
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

def define_gene_space(gene_list):
    gene_space = set()
    # Sort list of genes by start position
    gene_list.sort(key=lambda tup: tup[1])
    base_gene_list_position = 0
    comparison_gene_list_position = 1
    gene_of_interest = base_gene
    for i, feature in enumerate(gene_list[comparison_list_pos:]):
        base_gene = gene_list[base_gene_list_position]
        comparison_gene = gene_list[comparison_gene_list_position]
        base_gene_left = base_gene.location.start
        base_gene_right = base_gene.location.end
        comparison_gene_left = comparison_gene.location.start
        comparison_gene_right = comparison_gene.location.end
        base_gene_name = base_gene.feature.qualifiers['locus_tag'][0]
        comparison_gene_name = gene_list[comparison_gene_list_position].feature.qualifiers['locus_tag'][0]
        # If there's any gap between the two genes, add the gene as-is to the gene list
        if comparison_gene_left >= base_gene_left:
            gene_space.add(gene_list[base_gene_list_position])
        # Else if the comparison gene lies completely within the boundaries of the base gene,
        ### add the gene as-is to the gene list
        elif comparison_gene_left <= base_gene_right and comparison_gene_right <= base_gene_right:
            gene_space.add(gene_list[base_gene_list_position])
        # Else if the genes overlap, make a new compos
        else:

with open('229.opr', 'r') as operon_input:
    operon_reader = csv.reader(operon_input, delimiter = '\t')
    next(operon_reader, None) # ignore the header row
    input_list = []
    for row in operon_reader:
        input_list.append((row[0], row[2]))

record = SeqIO.read('NC_004350_2.gb', 'genbank')
reference_sequence = record.seq
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

# Loop through the operon dictionary and create a feature for each operon, then add these operons to a set.
operon_set = set()
for operon_id, gene_feature_list in operon_dict.items():
    start = gene_feature_list[0].location.start
    end = gene_feature_list[-1].location.end
    strand = gene_feature_list[0].strand
    new_operon_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type='operon')
    new_operon_feature.qualifiers['locus_tag'] = operon_id
    operon_set.add(new_operon_feature)
print('number of operons created: ', len(operon_set))

# Identify all the gene features that were in the original genbank but are NOT among those added to operons
gene_features_not_in_operons = complete_feature_set.difference(genes_in_operons)
print('gene features NOT in operons: ', len(gene_features_not_in_operons))

################################### DEFINE REGION 1: IGRs within operons ###################################

# Generate the gene list needed to feed the ss_IGR finder function (one list per operon)
igrs_within_operons = set()
for operon_id, gene_feature_list in operon_dict.items():
    operon_gene_list = feature_to_gene_list(gene_feature_list)
    # Search for IGRs within each operon.
    operon_igrs = get_ss_IGRs(reference_sequence, operon_gene_list, 1, 'intra_operon_IGR', strand, leading_trailing=False)
    igrs_within_operons.update(operon_igrs)
print('intra-operon IGRs identified: ', len(igrs_within_operons))

############################# DEFINE REGION 2: IGRs between transciptional units ############################

# Create a gene set consisting of operons plus genes *not* included in operons
transcriptional_units = operon_set.union(gene_features_not_in_operons)
print('number of transcriptional units: ', len(transcriptional_units))

# Separate transcriptional units on the forward and reverse strands
forward_tus = set()
reverse_tus = set()
for tu in transcriptional_units:
    if tu.strand == -1:
        reverse_tus.add(tu)
    elif tu.strand == 1:
        forward_tus.add(tu)
    else: print('STRAND ERROR: input features must be stranded')
print('number of forward TUs: ', len(forward_tus))
print('number of reverse TUs: ', len(reverse_tus))

# Create the list needed to feed the IGR finder
forward_tus_list = feature_to_gene_list(forward_tus)
reverse_tus_list = feature_to_gene_list(reverse_tus)

# Feed to IGR finder
forward_tu_igrs = get_ss_IGRs(reference_sequence, forward_tus_list, 0, 'inter_TU_IGR', 1, leading_trailing=True)
reverse_tu_igrs = get_ss_IGRs(reference_sequence, reverse_tus_list, 0, 'inter_TU_IGR', -1, leading_trailing=True)

# Reunite the sets into a single set of IGRs between transcriptional units.
igr_between_tus = forward_tu_igrs.union(reverse_tu_igrs)
print('number of IGRs in between TUs: ', len(igr_between_tus))

############################# DEFINE REGION 3: Coding regions (proteins and RNA: all genes) ############################

# Make a list of genes on forward and reverse strands
fwd_gene_list = []
rev_gene_list = []
for feature in record.features:
    if feature.type == 'gene':
        if strand == -1:
            rev_gene_list.append(feature)
        elif strand == 1:
            fwd_gene_list.append(feature)
        else: print('ERROR: not all genbank gene features stranded')


