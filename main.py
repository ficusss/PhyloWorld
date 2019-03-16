from phylo import PhyloManager

FILENAME_GENES = 'genes.txt'
FILENAME_SPECIES = 'cheburashkas.txt'
SEQ_FLODER = 'sequences'
ALIGNMENTS_FLODER = 'alignments'

def main():
    with open(FILENAME_GENES, 'r') as f:
        genes = f.read().split('\n')[:-1]
    with open(FILENAME_SPECIES, 'r') as f:
        species = f.read().split('\n')[:-1]
    
    pm = PhyloManager(genes, species, 'Homo sapien')
    pm.retrieve_sequences(SEQ_FLODER)
    #pm.alignments(SEQ_FLODER, ALIGNMENTS_FLODER)


if __name__ == '__main__':
    main()