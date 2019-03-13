from phyloutl import retrieve_sequences, merge_fasta, multiple_alignment
from phyloutl import find_dependent_rows, create_tree_from_seqs, create_tree_of_trees, create_tanglegram
from itertools import combinations
from ete3 import TreeStyle
from Bio import AlignIO
import glob
import os

variant = '4.2'
name = 'Mesheryakov Georgy'
email = 'herrberg@yandex.ru'
filename = 'Genes_2018.txt'
log_filename = 'logfile.txt'
fasta_folder = 'fasta'
dep_folder = 'dep'
phylo_folder = 'phylo_nj'
phylo_gr_folder = 'phylo_outgroup'
tangle_folder = 'tanglegrams'
species = """Pan troglodytes
Pan paniscus
Macaca mulatta
Mus musculus
Rattus norvegicus
Canis familiaris
Felis catus
Bos taurus
Gallus gallus
""".split('\n')[:-1]
outgroup_specie = 'Gallus gallus'.replace(' ', '_')
outgroup_gene = ''
families = [["Pan troglodytes", "Pan paniscus"],
            ['Rattus norvegicus', 'Mus musculus'],
            ["Canis familiaris", "Felis catus"]] # Expected families

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False

logfile = open(log_filename, 'w')
logfile.write("Student: {}, var. {}\n".format(name, variant))


with open(filename, 'r') as f:
    t = (line.split() for line in f.readlines())
    genes = [line[1] for line in t if line[-1] == variant]
    logfile.write("Genes: {}\n".format(' '.join(genes)))

# Alignment
d = retrieve_sequences(genes, species, logfile=logfile, email=email,
                       expected_families=families)
d_m = merge_fasta(d, logfile)
if not os.path.isdir('originals'):
    os.mkdir('originals')
for g in d_m:
    with open('{}/{}.fasta'.format('originals', g), 'w') as f:
        f.write(d_m[g])
logfile.close()


logfile = open(log_filename, 'w')
d_m = dict()
for g in genes:
    with open('{}/{}.fasta'.format('originals', g), 'r') as f:
        d_m[g] = f.read()
for gene, fasta in d_m.items():
    logfile.write("Running multiple alignment for gene {}...\n".format(gene))
    fasta_filename = '{}/{}.fasta'.format(fasta_folder, gene)
    algn = multiple_alignment(fasta, output_filename=fasta_filename)
logfile.close()


logfile = open(log_filename, 'w')
# Finding dependent rows
d_dep = dict()
if not os.path.isdir(dep_folder):
    os.mkdir(dep_folder)
for gene in genes:
    aln = AlignIO.read("fasta/{}.fasta".format(gene), "fasta")
    d_dep[gene] = find_dependent_rows(aln, equals=True)
    with open('dep/{}.txt'.format(gene), 'w') as f:
        for key, t in d_dep[gene].items():
            f.write('{}: ({})\n{}\n'.format(key, len(t), ','.join((str(i) for i in t))))
logfile.close()


logfile = open(log_filename, 'w')
# Creating phylogenetic trees
for gene in genes:
    logfile.write("Creating a phylogenetic tree for gene {}...\n".format(gene))
    phylo_filename = "{}/{}.nwk".format(phylo_folder, gene)
    fasta_filename = "{}/{}.fasta".format(fasta_folder, gene)
    tree = create_tree_from_seqs(fasta_filename, phylo_filename,
                                 outgroup=outgroup_specie)
    tree.render("{}/{}.png".format(phylo_folder, gene), tree_style=ts, dpi=600)
logfile.close()

logfile = open(log_filename, 'w')
# Creating phylogenetic trees with outgroup present
for gene in genes:
    logfile.write("Creating a phylogenetic tree for gene {}...\n".format(gene))
    phylo_filename = "{}/{}.nwk".format(phylo_gr_folder, gene)
    fasta_filename = "{}/{}.fasta".format(fasta_folder, gene)
    tree = create_tree_from_seqs(fasta_filename, phylo_filename,
                                 outgroup=outgroup_specie,
                                 delete_outgroup=False)
    tree.render("{}/{}.png".format(phylo_folder, gene), tree_style=ts, dpi=600)
logfile.close()


logfile = open(log_filename, 'w')
# Creating a tree of trees
logfile.write("Creating a tree of trees...\n")
trees = {filename.split('/')[-1].split('.')[-2]: filename
         for filename in glob.glob('{}/*.nwk'.format(phylo_folder))} 
tree_of_trees = create_tree_of_trees(trees, '{}/tot.nwk'.format(phylo_folder),
                                     outgroup=outgroup_gene)
tree_of_trees.render("{}/tot.png".format(phylo_folder), tree_style=ts, dpi=600)
trees = {filename.split('/')[-1].split('.')[-2]: filename
         for filename in glob.glob('{}/*.nwk'.format(phylo_gr_folder))} 
tree_of_trees = create_tree_of_trees(trees,
                                     '{}/tot.nwk'.format(phylo_gr_folder),
                                     outgroup=outgroup_gene)
tree_of_trees.render("{}/tot.png".format(phylo_gr_folder), tree_style=ts,
                     dpi=600)
logfile.close()


logfile = open(log_filename, 'w')   
# Creating dendrograms
logfile.write("Creating dendrograms via R script...\n")
for a, b in combinations(genes, 2):
    filename_a = '{}/{}.nwk'.format(phylo_folder, a)
    filename_b = '{}/{}.nwk'.format(phylo_folder, b)
    out = '{}/{}_{}.png'.format(tangle_folder, a, b)
    create_tanglegram(filename_a, filename_b, out, a, b)
logfile.close()
