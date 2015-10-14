#from ete2 import NCBITaxa
ncbi = NCBITaxa()

def get_closest(taxids, target):
    ''' Get the closest taxonomy ID to the target ID based on tree distance.

    Args:
        taxids (iterable): An iterable of taxonomy IDs.
        target (int): The taxonomy ID to use a reference.

    Returns:
        (int, int): The distance and closest taxonomy ID.

    '''

    dist = None
    closest = None

    L_T = ncbi.get_lineage(target)[::-1]
    #print 'Lineage for {0}: {1}'.format(target, L_T)

    for S_i in taxids:
        L_si = ncbi.get_lineage(S_i)[::-1]
        #print 'Checking {0} with lineage {1}'.format(S_i, L_si)
        for n, s in enumerate(L_si):
            if s in L_T:
                d = n + L_T.index(s)
                if dist is None or d < dist:
                    #print 'Found new distance {0} < {1}'.format(d, dist)
                    dist = d
                    closest = S_i
                    break
    return dist, closest

def get_closest_per_ortho(df, target):
    ''' Find the best hit out of the closest hit for each query.

    Args:
        df (DataFrame): The hits, in MAF format.
        target (int): The taxonomy ID to use as a reference.

    Returns:
        (DataFrame): The filtered hits.
    '''

    def best(group):
        dist, closest = get_closest(group.odb8_level, target)
        closest_subset = group.query('odb8_level == {0}'.format(closest))
        closest_subset.sort('bitscore', inplace=True)
        try:
            return closest_subset.iloc[0]
        except IndexError as e:
            print group
            print dist, closest
            raise

    top_per_prot = df.groupby(['q_name', 'odb8_prot_id']).apply(best)
    top_per_query = df.groupby('q_name').apply(best)

    return top_per_query.reset_index()


