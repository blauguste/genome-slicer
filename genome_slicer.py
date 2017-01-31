from Bio import SeqIO, SeqRecord, Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
from collections import defaultdict
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from full_genome_promoters import find_35_10_promoters

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

def union(feature_list, strand):
    gene_coverage_feature_set = set()
    intervals = []
    for feature in feature_list:
        intervals.append((feature.location.start, feature.location.end))
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
        # test for intersection between lower and higher:
        # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    for entry in merged:
        new_feature = make_IGR_feature(entry[0], entry[1], 'gene_no_overlap', strand)
        gene_coverage_feature_set.add(new_feature)
    return gene_coverage_feature_set

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
    strand = operon_gene_list[0][3]
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
        if feature.strand == -1:
            rev_gene_list.append(feature)
        elif feature.strand == 1:
            fwd_gene_list.append(feature)
        else: print('ERROR: not all genbank gene features stranded')

# Find the union between all the location sets
fwd_coverage_set = union(fwd_gene_list, 1)
rev_coverage_set = union(rev_gene_list, -1)
print('forward features reduced from %s genes to %s coverage areas' % (str(len(fwd_gene_list)), str(len(fwd_coverage_set))))
print('reverse features reduced from %s genes to %s coverage areas' % (str(len(rev_gene_list)), str(len(rev_coverage_set))))

# Reunite the stranded coverage sets into a single set (both strands)
complete_coverage_set = fwd_coverage_set.union(rev_coverage_set)
print('length of complete coverage set: ', len(complete_coverage_set))
print('bp in complete coverage set: ', sum(len(feature) for feature in complete_coverage_set))

############################# DEFINE REGION 4: Regions opposite coding regions ############################

# Create an identical set of features to complete coverage feature, but on the opposite strand
complementary_coverage_set = set() 
for feature in complete_coverage_set:
    start = feature.location.start
    end = feature.location.end
    complementary_strand = feature.strand * -1
    complementary_location = FeatureLocation(start, end, complementary_strand)
    complementary_feature = SeqFeature(complementary_location, type='opposite_gene_region', strand=complementary_strand)
    complementary_coverage_set.add(complementary_feature)

print('length of complementary coverage set: ', len(complementary_coverage_set))

###########################################################################################################

#### need to prioritize which categories take precedence, and remove the overlaps!
# 1. coding regions (region 3)
# 2. IGRs within operons (region 1)
# 3. Area opposite genes (region 4)
# 4. IGRs between transcriptional units (region 2)

bp_region1 = sum(len(feature) for feature in igrs_within_operons) # IGRs within operons
bp_region2 = sum(len(feature) for feature in igr_between_tus) # IGRs between transcriptional units
bp_region3 = sum(len(feature) for feature in complete_coverage_set) # Gene coverage area (protein-coding & RNA)
bp_region4 = sum(len(feature) for feature in complementary_coverage_set) # Area opposite genes (region 3)

total_bp = bp_region1 + bp_region2 + bp_region3 + bp_region4
print('total base pairs: ', total_bp)

print(bp_region1, bp_region2, bp_region3, bp_region4)

bp_fwd_tus = sum(len(feature) for feature in forward_tus) # Transcriptional units on the forward strand
bp_rev_tus = sum(len(feature) for feature in reverse_tus) # Transcriptional units on the reverse strand
print('bp forward tus: ', bp_fwd_tus)
print('bp reverse tus: ', bp_rev_tus)

sum_len_gene_features = sum(len(feature) for feature in record.features)
print('sum of all gene feature length: ', sum_len_gene_features)
sum_fwd_coverage = sum(len(feature) for feature in fwd_coverage_set)
print('sum of fwd gene coverage: ', sum_fwd_coverage)
sum_rev_coverage = sum(len(feature) for feature in rev_coverage_set)
print('sum of fwd gene coverage: ', sum_rev_coverage)

################### CHECK WHICH SETS PROMOTER POSITIONS ARE IN (i.e. where are the promoters?) ###############

# Create ranges tuples to search for each set, forward and reverse
#forward_dict = {}
#forward_dict['igr_within_operon'] = list((feature.location.start.position, feature.location.end.position) for feature in igrs_within_operons if feature.strand == 1)
#forward_dict['igr_between_tu'] = list((feature.location.start.position, feature.location.end.position) for feature in igr_between_tus if feature.strand == 1)
#forward_dict['gene_coding_region'] = list((feature.location.start.position, feature.location.end.position) for feature in complete_coverage_set if feature.strand == 1)
#forward_dict['opposite_gene_coding'] = list((feature.location.start.position, feature.location.end.position) for feature in complementary_coverage_set if feature.strand == 1)

#reverse_dict = {}
#reverse_dict['igr_within_operon'] = list((feature.location.start.position, feature.location.end.position) for feature in igrs_within_operons if feature.strand == -1)
#reverse_dict['igr_between_tu'] = list((feature.location.start.position, feature.location.end.position) for feature in igr_between_tus if feature.strand == -1)
#reverse_dict['gene_coding_region'] = list((feature.location.start.position, feature.location.end.position) for feature in complete_coverage_set if feature.strand == -1)
#reverse_dict['opposite_gene_coding'] = list((feature.location.start.position, feature.location.end.position) for feature in complementary_coverage_set if feature.strand == -1)

promoters_of_interest = find_35_10_promoters('NC_004350_2.gb', 'strep_mutans_ext_35_only_test.csv')

with open('TSS_regions_test.csv', 'w') as TSS_output:
    writer = csv.writer(TSS_output, delimiter = '\t')
    for promoter in promoters_of_interest:
        matching_regions = []
        promoter_pos = promoter[1]
        if promoter[0] == 'F':
            strand = 1
            TSS_pos = promoter_pos + 13
        if promoter[0] == 'R':
            strand = -1 
            TSS_pos = promoter_pos - 13
        for feature in igrs_within_operons:
            if feature.strand == strand:
                if TSS_pos in feature:
                    matching_regions.append('IGRs_within_operons')
        for feature in igr_between_tus:
            if feature.strand == strand:
                if TSS_pos in feature:
                    matching_regions.append('IGRs_between_TUs')
        for feature in complete_coverage_set:
            if feature.strand == strand:
                if TSS_pos in feature:
                    matching_regions.append('coding_region')
        for feature in complementary_coverage_set:
            if feature.strand == strand:
                if TSS_pos in feature:
                    matching_regions.append('complementary_to_coding_region')
        writer.writerow([strand, promoter_pos, TSS_pos, matching_regions])

############################ VISUALIZATION!! #######################################################

for new_feature in igrs_within_operons:
    record.features.append(new_feature) # type: #intra_operon_IGR
for new_feature in igr_between_tus:
    record.features.append(new_feature) # type: #inter_TU_IGR
for new_feature in complete_coverage_set:
    record.features.append(new_feature) # type: #gene_no_overlap

outpath = 'strep_mutans_sliced.gb'
SeqIO.write(record, open(outpath, 'w'), 'genbank')

gd_diagram = GenomeDiagram.Diagram('Strep mutans, sliced!')
gd_track_for_features = gd_diagram.new_track(1, name='newly sliced features')
gd_feature_set = gd_track_for_features.new_set()

for gb_feature in record.features:
    if gb_feature.type == 'intra_operon_IGR':
        color = colors.orangered
        gd_feature_set.add_feature(gb_feature, color=color)
    if gb_feature.type == 'inter_TU_IGR':
        color = colors.darkslateblue
        gd_feature_set.add_feature(gb_feature, color=color)
    if gb_feature.type == 'gene_no_overlap':
        color = colors.limegreen
        gd_feature_set.add_feature(gb_feature, color=color)
    
gd_diagram.draw(format='circular', circular=True, pagesize=(20*cm,20*cm), start=0, end=len(record), circle_core=0.7)
gd_diagram.write('strep_mutans_sliced.pdf', 'PDF')
