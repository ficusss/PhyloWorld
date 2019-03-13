import sys, os, glob
import numpy as np
import subprocess
import ete3
from Bio import Entrez
from time import sleep
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Align import MultipleSeqAlignment as MSA
import Bio.Blast.NCBIXML as BlastXML
import Bio.Blast.NCBIWWW as Blast
from Bio.Seq import Seq
from skbio import nj, DistanceMatrix
from Bio import AlignIO, SeqIO
from itertools import permutations
from collections import defaultdict, deque
from itertools import combinations, product
from pandas import DataFrame

__rscript_tangle = '''library(phytools)
library(dendextend)
library(ape)
args = commandArgs(trailingOnly=TRUE)
a <- args[1]
b <- args[2]
a <- chronos(read.newick(a))
b <- chronos(read.newick(b))

png(filename=args[3], width=768)
tanglegram(a, b, main_left=args[4], main_right=args[5], axes=FALSE, sort=TRUE, margin_outer=10, margin_inner=8)
dev.off()'''

__rscript_distance = '''library(phytools)
library(ape)
args = commandArgs(trailingOnly=TRUE)
a <- args[1]
b <- args[2]
c <- args[3]
a <- read.newick(a)
b <- read.newick(b)
print(dist.topo(a, b, method=c))
'''

__RETRIES = 3
__SLEEP_TIME_RETRY = 2


def parse_path(out):
    t = out.split('/')
    output_dir = str()
    if len(t) == 1:
        output_filename = t[0]
    else:
        output_filename = t[-1]
        for folder in t[:-1]:
            output_dir += folder
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            output_dir += '/'
    return output_dir, output_filename


def retry(f, logfile, *args, **argw):
    global __RETRIES
    global __SLEEP_TIME_RETRY
    last_error = None
    for i in range(__RETRIES):
        try:
            return f(*args, **argw)
        except Exception as e:
            if logfile:
                logfile.write("Error: {}. Attempt {}/{}\n".format(e, i + 1,
                              __RETRIES))
            last_error = e
            sleep(__SLEEP_TIME_RETRY)
    raise last_error

def multiple_alignment(fastas, output_filename='align.fasta',
                       tmp_filename='tmp.fasta', codon_preservation=True,
                       **args):
    def run_aln(in_file, out_file, **args):
        ClustalwCommandline('clustalw2', infile=in_file,
                        outfile=out_file + '.aln', **args)()
        align = AlignIO.read(out_file + '.aln', 'clustal')
        os.remove(out_file + '.aln')
        with open(out_file, "w") as f:
            AlignIO.write(align, f, 'fasta')
        return align
    parse_path(output_filename)
    if type(fastas) is not str:
        fastas = '\n'.join(fastas)
    if '>' in fastas:
        with open(tmp_filename, 'w') as f:
            f.write(fastas)
        fastas = tmp_filename
        tmp_filename = True
    else:
        tmp_filename = False
    if not codon_preservation:
        align = run_aln(fastas, output_filename)
    else:
        with open(fastas, 'r') as f:
            seqs = SeqIO.parse(f, 'fasta')
            d = {seq.name: seq for seq in seqs}
        tmp = fastas + '.amino.tmp'
        s = '\n'.join('>{}\n{}'.format(name, str(d[name].translate(cds=True).seq))
                      for name in d)
        with open(tmp, 'w') as f:
            f.write(s)
        tmp_out = fastas + '.amino.fasta'
        alns = run_aln(tmp, tmp_out)
        os.remove(tmp)
        os.remove(tmp_out)
        d_aln = {aln.id: str(aln.seq) for aln in alns}
        out_text = list()
        for name, seq in d.items():
            t = str(seq.seq)
            codons = (t[i:i + 3] for i in range(0, len(t), 3))
            reconstructed_nt = str()
            aminos = d_aln[name]
            for amino in aminos:
                if amino == '-':
                    reconstructed_nt += '-' * 3
                else:
                    reconstructed_nt += next(codons)
            reconstructed_nt += next(codons)
            out_text.append('>{}\n{}'.format(name, reconstructed_nt))
        with open(output_filename, 'w') as f:
            f.write('\n'.join(out_text))
        align = AlignIO.read(output_filename, 'fasta')
    if tmp_filename:
        os.remove(fastas)
    return align


def find_dependent_rows(alignments, equals=True):
    def check(x, y):
        src, target = list(), list()
        for a, b in zip(x, y):
            try:
                i = src.index(a)
                if b != target[i]:
                    return False
            except ValueError:
                if b in target:
                    return False
                src.append(a)
                target.append(b)
        return True
    if type(alignments) is MSA:
        seqs = np.array([list(s.seq) for s in alignments]).T
    else:
        seqs = np.array([list(s) for s in alignments]).T
    d = dict()
    seqs = list(np.apply_along_axis(lambda r: ''.join(r), axis=1, arr=seqs))
    uniqs = list(len({letter for letter in col}) if '-' not in col else 0
                 for col in seqs)
    visited = list(False for _ in range(len(seqs)))
    m = len(seqs[0])
    groups = [deque() for _ in range(m)]
    for i, k in enumerate(uniqs):
        if k:
            groups[k - 1].append(i)
    d = dict()
    for i, (col, uniq) in enumerate(zip(seqs, uniqs)):
        if uniq and not visited[i]:
            d[col] = list()
            group = groups[uniq - 1]
            for j in group:
                if not visited[j] and check(col, seqs[j]):
                    if equals or col != seqs[j]:
                        d[col].append(j)
                    visited[j] = True
        visited[i] = True
        for group in groups:
            while group and (visited[group[0]] or group[0] <= i):
                group.popleft()
    return {key: items for key, items in d.items() if len(items) > 1}
       

def create_tree_from_seqs(filename, out, outgroup='', model='GTRGAMMA',
                          tail_deletion=False, n=40, delete_outgroup=True):
    output_dir, output_filename = parse_path(out)
    if tail_deletion:
        regions = list()
        alns = AlignIO.read(filename, 'fasta')
        a = None
        for i in range(len(alns[0])):
            complete = True
            for aln in alns:
                if aln[i] == '-':
                    complete = False
                    break
            if complete:
                if a is None:
                    a = i
            else:
                if a is not None:
                    regions.append((a, i))
                a = None
        if a is not None:
            regions.append((a, i))
        s, _ = max(regions, key=lambda reg: (reg[1] - reg[0]) / (reg[0] + 1))
        alns = MSA([aln[s:] for aln in alns])
        filename = filename + '.tmp.fasta'
        AlignIO.write(alns, filename, 'fasta')           
    RaxmlCommandline(sequences=filename, model=model,
                     parsimony_seed=2018, num_replicates=n,
                     name=output_filename)()
    if tail_deletion:
        os.remove(filename)
    for filename in glob.glob("RAxML_*"):
        if output_filename in filename:
            if not filename.endswith('bestTree.' + output_filename):
                os.remove(filename)
            else:
                os.rename(filename, out)
    if os.path.isfile('tmp.dnd'):
        os.remove('tmp.dnd')
    tree = ete3.Tree(out)
    os.remove(out)
    if outgroup is not None:
        outgroup_tree(tree, outgroup, delete_outgroup)
    tree.write(outfile=out)
    return ete3.Tree(out)

def get_dist_trees(ref: str, targets: list, metric='score'):
    if metric == 'RF' or metric == 'nRF':
        result = subprocess.check_output(['ete3', 'compare', '-t'] +
                                         targets + ['-r', ref,
                                                    '--taboutput']).decode('UTF-8')
        result = result.rstrip('\n')
        lines = [t.split('\t') for t in result.split('\n')]
        i = lines[0].index(metric)
        return [float(line[i]) for line in lines[1:]]
    else:
        global __rscript_distance
        tmp_filename = 'tmp-filename.r'
        with open(tmp_filename, 'w') as f:
            f.write(__rscript_distance)
        def get_dist(a, b, metric):
            res = subprocess.check_output(['Rscript', tmp_filename, a, b,
                                           metric]).decode('UTF-8').split('\n')
            res = [t for t in res if len(t)]
            res = res[-1].split(' ')[-1]
            print(res)
            return float(res)
        ret = [get_dist(ref, b, metric) for b in targets]
        os.remove(tmp_filename)
        return ret

def create_tanglegram(filename_a, filename_b, output, title_a, title_b):
    global __rscript_tangle
    parse_path(output)
    tmp_filename = 'tmp-tangle.r'
    with open(tmp_filename, 'w') as f:
        f.write(__rscript_tangle)
    subprocess.check_output(['Rscript', tmp_filename,
                             filename_a, filename_b, output,
                             title_a, title_b])
    os.remove(tmp_filename)

def outgroup_tree(tree, outgroup='', delete_outgroup=False):
    if outgroup:
        for leaf in tree.get_leaves():
            if leaf.name == outgroup:
                outgroup = leaf
                break
    else:
        outgroup = tree.get_midpoint_outgroup()
    tree.set_outgroup(outgroup)
    if delete_outgroup:
        outgroup.delete()

def create_tree_of_trees(trees: dict, out, metric='score', outgroup=''):
    output_dir, output_filename = parse_path(out)
    names = sorted(trees.keys())
    n = len(names)
    m = np.zeros((n, n))
    for i, ref in enumerate(names[:-1]):
        targets = [trees[names[j]] for j in range(i + 1, n)]
        dists = get_dist_trees(trees[ref], targets, metric)
        m[i + 1: n, i] = dists
    m += m.T
    dist_m = DistanceMatrix(m, names)
    tree = ete3.Tree(nj(dist_m, result_constructor=str))
    if outgroup is not None:
        outgroup_tree(tree)
    if out:
        tree.write(outfile=out)
    return tree

def get_gene(gene, specie, logfile=None, db='nuccore'):
    search_str = "{}[Gene] AND {}[Organism]".format(gene, specie)
    search_str += " AND mRNA[Filter] AND RefSeq[Filter]"
    h = retry(Entrez.esearch, logfile, term=search_str, db=db)
    records = Entrez.read(h)
    fasta, i = None, None
    for i in records['IdList']:
        fasta = get_gene_by_id(i, logfile, db=db)
        if fasta:
            break
    if i is None:
        s = "Incorrect command for Entrez.esearch: {}".format(search_str)
        raise Exception(s)
    return i, fasta

def get_gene_by_id(identificator, logfile=None, db='nuccore'):
    h = retry(Entrez.efetch, logfile, db='nuccore', id=identificator,
              rettype="fasta_cds_na", retmode="text")
    record = h.read()
    fasta = record.rstrip('\n')
    t = fasta.split('\n')
    if len(t) < 2 or len(t[1].replace(' ', '')) <= 1:
        if logfile:
            s = "Couldn't retrieve a CDS from {}.\n".format(identificator)
            logfile.write(s)
        return None
    return fasta

def get_homologous_seqs(nuc_fasta, specie, logfile=None, ret_max=5):
    t = nuc_fasta.split('\n')
    h = t[0]
    t = ''.join(t[1:])
    amino_seq = str(Seq(t).translate(cds=True))
    amino_fasta = '{}\n{}'.format(h, amino_seq)
    blast_res = retry(Blast.qblast, logfile, 'tblastn', 'nt', amino_fasta,
                      expect=0.1, entrez_query="{}[ORGANISM]".format(specie))
    rec = BlastXML.read(blast_res)
    if not rec.alignments:
        if logfile:
            s = "tblastn found no results, trying blastn (no megablast).\n"
            logfile.write(s)
        blast_res = retry(Blast.qblast, logfile, 'blastn', 'nt', nuc_fasta,
                          expect=0.1,
                          entrez_query="{}[ORGANISM]".format(specie))
        rec = BlastXML.read(blast_res)
    fastas = list()
    for aln in rec.alignments:
        i = int(aln.title.split('|')[1])
        if logfile:
            s = "Sending a request for gene ID = {}...\n".format(i)
            logfile.write(s)
        fasta = get_gene_by_id(i, logfile)
        if fasta:
            fastas.append(fasta)
            if logfile:
                logfile.write("Success.\n")
            if len(fastas) >= ret_max:
                break
    return fastas

def get_homologous_seqs_retry(nuc_fasta, specie, logfile=None, ret_max=5):
    def wrapper():
        ret = get_homologous_seqs(nuc_fasta, specie, logfile, ret_max)
        if not ret:
            raise Exception("Couldn't find any homologic seqs for {}".format(specie))
        return ret
    return retry(wrapper, logfile)

def get_short_headed(fasta, specie):
    header = '>{}'.format(specie.replace(' ', '_'))
    fasta = '\n'.join([header] + fasta.split('\n')[1:])
    return fasta

def naive_score(seq_a, seq_b):
    return sum(1 for a, b in zip(seq_a, seq_b) if a == b)    

def retrieve_sequences(genes, species, ref_organism='Homo sapien',
                       logfile=sys.stdout, email='herrberg@yandex.ru',
                       short_headers=True, expected_families=list()):
    Entrez.email = email
    d = dict()
    for gene in genes:
        print("{}".format(gene))
        logfile.write("Sending request to NCBI for  gene {}...\n".format(gene))
        i, fasta = get_gene(gene, ref_organism, logfile)
        logfile.write("ID = {}\n".format(i))
        if short_headers:
                    fasta = get_short_headed(fasta, ref_organism)
        d[gene] = [fasta]
        tmp = dict()
        for specie in species:
            print("\t{}".format(specie))
            logfile.write("Sending a BLAST request to NCBI for specie {}\n".format(specie))
            res = get_homologous_seqs_retry(d[gene][0], specie, logfile)
            if not res:
                s = "Couldn't find anything for specie {}.\n".format(specie)
                logfile.write(s)
            else:
                tmp[specie] = res
        # Family testing
        logfile.write("Performing family testing...")
        tmp_seq = dict()
        for specie, fastas in tmp.items():
            seqs = [''.join(fasta.split('\n')[1:]) for fasta in fastas]
            tmp_seq[specie] = seqs
        for fam in expected_families:
            for specie_a, specie_b in combinations(fam, 2):
                if specie_a not in tmp or specie_b not in tmp:
                    print("Warning: couldn't retrieve some of the species.")
                    continue
                a = tmp_seq[specie_a]
                inds_a = list(range(len(a)))
                b = tmp_seq[specie_b]
                inds_b = list(range(len(b)))
                i, j = max(product(inds_a, inds_b),
                           key=lambda x: naive_score(a[x[0]], b[x[1]]))
                t = [tmp[specie_a][i]]
                tmp[specie_a] = t
                tmp_seq[specie_a] = t
                t = [tmp[specie_b][j]]
                tmp[specie_b] = t
                tmp_seq[specie_b] = t
        logfile.write('Selected IDs:\n')
        for specie, items in tmp.items():
            fasta = items[0]
            logfile.write('{}: {}\n'.format(specie, fasta.split('\n')[0][1:]))
            if short_headers:
                fasta = get_short_headed(fasta, specie)
            d[gene].append(fasta)
    return d    

def get_reverse_evol_nucs(tree: ete3.Tree):
    pass

def merge_fasta(dictionary, logfile=sys.stdout):
    d = dict()
    logfile.write("Merging FAST data into multi-FAST format...\n")
    for gene, fastas in dictionary.items():
        d[gene] = '\n'.join(fastas)
    return d