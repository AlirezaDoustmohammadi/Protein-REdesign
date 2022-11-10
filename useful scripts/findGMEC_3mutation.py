import osprey


def running_osprey(candidate_amino_acids, candidate_amino_acids_fixed_pos, protein, chain, mutation_required_pos,
                   mutation_required_aa, fixed_mutation_required_pos, n_core):

    osprey.start()

    parallelism = osprey.Parallelism(cpuCores=n_core)

    # calculates
    for i in range(0, len(mutation_required_aa)):
        #  Candidate amino acids for this position
        candidate_amino_acids_this_pos = candidate_amino_acids[:]
        candidate_amino_acids_this_pos.remove(mutation_required_aa[i])
        if mutation_required_aa[i] == 'ASP':
            candidate_amino_acids_this_pos.remove('PRO')
        strand = osprey.Strand(protein + ".pdb")
        # dynamic position
        strand.flexibility[chain + mutation_required_pos[i]].setLibraryRotamers(*candidate_amino_acids_this_pos)
        # fixed positions
        strand.flexibility[chain + fixed_mutation_required_pos[0]].setLibraryRotamers(*candidate_amino_acids_fixed_pos)
        strand.flexibility[chain + fixed_mutation_required_pos[1]].setLibraryRotamers(*candidate_amino_acids_fixed_pos)
        pdb_file = protein + '_' + chain + '_' + mutation_required_pos[i] + '.pdb'

        # make the conf space
        confSpace = osprey.ConfSpace(strand)

        # choose a forcefield
        ffparams = osprey.ForcefieldParams()

        # how should we compute energies of molecules?
        ecalc = osprey.EnergyCalculator(confSpace, ffparams, parallelism=parallelism)

        # how should we define energies of conformations?
        confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

        # how should confs be ordered and searched?
        emat = osprey.EnergyMatrix(confEcalc)
        astar = osprey.AStarMPLP(emat, confSpace)

        # find the best sequence and rotamers
        gmec = osprey.GMECFinder(astar, confEcalc).find()

        # write the rigid GMEC to a pdb
        gmecStructure = confSpace.makeMolecule(gmec.getAssignments())
        osprey.writePdb(gmecStructure, pdb_file)


if __name__ == '__main__':

    # number of cpu core
    n_core = 9

    # Candidate amino acids that can be substituted for natural amino acids in proteins
    # for variable positions
    candidate_amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                             'LEU', 'MET', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'PRO']

    # for fix position
    candidate_amino_acids_fixed_pos = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                                       'LEU', 'MET', 'ASN', 'GLN', 'ARG', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

    # define protein
    protein = '3hqr'
    chain = 'S'
    # Mutation-required positions
    mutation_required_pos = ['558', '559', '560', '561', '562', '563', '565', '566', '568',
                             '569', '570', '571', '572', '573', '574']
    # The amino acids that corresponds to the positions
    mutation_required_aa = ['ASP', 'LEU', 'GLU', 'MET', 'LEU', 'ALA', 'TYR', 'ILE', 'MET', 'ASP', 'ASP',
                            'ASP', 'PHE', 'GLN', 'LEU']

    # fixed mutation-required positions
    fixed_mutation_required_pos = ['564', '567']

    running_osprey(candidate_amino_acids, candidate_amino_acids_fixed_pos, protein, chain, mutation_required_pos,
                   mutation_required_aa, fixed_mutation_required_pos)
