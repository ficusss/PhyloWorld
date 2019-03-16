import multiprocessing as mp
import phyloutl as phy
import os


def get_gene(organism: str, gene: str):
    return phy.get_short_headed(phy.get_gene(gene, organism)[1], organism)

class PhyloWorker(object):
    """The worker retrieves homologous seqs from BLAST.
    Keyword arguments:
        gene -- The name of gene.
        reference -- Either a name of referense organism i.e. 'Homo sapien',
                     or a FASTA sequence.
    """
    def __init__(self, gene: str, reference):
        if reference.split('\n')[0].strip(' ')[0] != '>':
            reference = get_gene(reference, gene)
        self.ref = reference
        self.gene = gene

    def get_homologs(self, specie, ret_max=4):
        """Retrieve homologs for specie.
        Keyword arguments:
            specie  -- The specie's name.
            ret_max -- Maximum number of homologs to return.
        Returns:
            A list of homologous genes' FASTA CDS.
        """
        return phy.get_homologous_seqs(self.ref, specie, ret_max=ret_max)


class PhyloManager(object):
    """Multithreads a set of phylo_workes each worker per gene.
    """
    def __init__(self, genes, species, ref_organism='Homo sapien'):
        self.ref = ref_organism
        self.genes = genes
        self.species = species

    def retrieve_sequences(self, folder: str):
        """Retrieves sequences from GenBank/BLAST db and saves them to folder.
        Keyword arguments:
            folder            -- Folder name.
        """
        if not os.path.isdir(folder):
            os.mkdir(folder)
        processes = [mp.Process(target=self.run_seqs_retrievement,
                                args=(gene, folder))
                     for gene in self.genes]
        for proc in processes:
            proc.start()

    def alignments(self, infolder: str, outfolder: str):
        """Performs a MSA on genes present in infolder and saves results to
           outfolder.
        Keyword argument:
            infolder  -- The name of folder where .fasta files are stored.
            outfolder -- The name of folder where alignments are to be stored.
        """
        if not os.path.isdir(outfolder):
            os.mkdir(outfolder)
        processes = [mp.Process(target=self.run_seqs_alignment,
                                args=(gene, infolder, outfolder))
                     for gene in self.genes]
        for proc in processes:
            proc.start()

    def run_seqs_retrievement(self, gene, folder):
        ext = '{}/{}-extra'.format(folder, gene)
        if not os.path.isdir(ext):
            os.mkdir(ext)
        worker = PhyloWorker(gene, self.ref)
        d = {specie: worker.get_homologs(specie) for specie in self.species}
        fasta_lt = list()
        for specie, fastas in d.items():
            it = iter(fastas)
            fasta_lt.append(next(it))
            for i, fasta in enumerate(it):
                with open('{}/{}-{}.fasta'.format(ext, specie, i), 'w') as f:
                    f.write(fasta)
        merged = '\n'.join(fasta_lt)
        with open('{}/{}.fasta'.format(folder, gene), 'w') as f:
            f.write(merged)

    def run_seqs_alignment(self, gene, infolder, outfolder):
        with open('{}/{}.fasta'.format(infolder, gene), 'r') as f:
            fastas = f.read()
        out = '{}/{}.fasta'.format(outfolder, gene)
        out_tmp = '{}/{}.fasta.tmp'.format(outfolder, gene)
        phy.multiple_alignment(fastas, output_filename=out,
                               tmp_filename=out_tmp)
