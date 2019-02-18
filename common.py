

def check_species(species):
    species = sys.argv[1]
    if species not in ('mm10', 'hg38'):
        print('Species "%s" not found' % species)
        print('Valid Species codes are:')
        print('    hg38 - human')
        print('    mm10 - mouse')
        print()
        return False
    return True
