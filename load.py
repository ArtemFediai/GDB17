# load and get 100 mols (or any other other number like 1000) of every atoms count
import os.path
from decimal import Decimal
from rdkit import Chem
from rdkit import RDLogger
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn
import pickle
import os

RDLogger.DisableLog('rdApp.*')

DATA_FOLDER = '/home/artem/dataset_gdb17'
DATA_FILE_NAME = 'GDB17.50000000.smi'
DATA_PATH = DATA_FOLDER + '/' + DATA_FILE_NAME
# MOL_NUM_LIMIT = False  # if False, the whole data set
MOL_SIZES_OF_INTEREST = [10, 11, 12, 13, 14, 15, 16, 17, 18]  # number of atoms
SAMPLE_SIZE = 100000


def main():
    print('I run main')

    # atoms_count_list = []
    # unique_atom_counts = []
    # occurences_number = []
    # permuted_numbers = []

    saved_dict = {
        'atom_count_list.pkl': [],
        'unique_atom_counts.pkl': [],
        'occurences_number.pkl': [],
        'permuted_numbers.pkl': []
    }

    for save_keys in saved_dict:
        with open(DATA_FOLDER + '/' + save_keys, 'rb') as fid:
            saved_dict[save_keys] = pickle.load(fid)

    atoms_count_list = saved_dict['atom_count_list.pkl']  #
    unique_atom_counts = saved_dict['unique_atom_counts.pkl']  # number of atoms
    occurences_number = saved_dict['occurences_number.pkl']  # 0 atoms, 1 atom, 2 atoms, ...
    permuted_numbers = saved_dict['permuted_numbers.pkl']

    # show counts of molecules of every size
    ax = seaborn.countplot(atoms_count_list)
    plt.show()

    # print(f'Number of atoms per molecules in the dataset: {unique_atom_counts}')
    # print(f'Count of molecules with respecting atom number: {occurences_number}')

    # atoms_num: mol_count, i.e.: 10:3379, 11:48909, 12:948094, ...
    dict_atoms_num_mol_count = dict(zip(np.arange(len(occurences_number)), occurences_number))
    print(dict_atoms_num_mol_count)

    # ids of molecules of every size, dictionary
    # 10: [344, 1244, 334554, 2344, 3355, 33455, 34355, ...], 11: [433, 445, 3354566, ...], ...
    dict_num_atoms_mol_idx = {}
    for mol_size_iterate in MOL_SIZES_OF_INTEREST:
        print(f'mol size of interest = {mol_size_iterate}')
        dict_num_atoms_mol_idx[mol_size_iterate] = \
            [mol_idx for mol_idx, mol_size in enumerate(atoms_count_list) if mol_size == mol_size_iterate]

    current_path = DATA_FOLDER + '/' + 'dict_num_atoms_mol_idx.pkl'  # <-- main dict. Can be used in a separate file
    with open(current_path, 'wb') as fid:
        pickle.dump(dict_num_atoms_mol_idx, fid)

    # from this point on, it may be reasonable to make a new file. dict_num_atoms_mol_idx.pkl is mostly needed
    ####################################################################################################################
    # dict that contains mol idx of random molecules with X number of atoms. X = 10...18 (why 18 is included, no idea)
    dict_100 = {}
    for mol_size_iterate in MOL_SIZES_OF_INTEREST:
        print(f'mol_size by generating random molecules: {mol_size_iterate}')
        my_population = dict_num_atoms_mol_idx[mol_size_iterate]
        # if len(my_population) > SAMPLE_SIZE:
        try:
            dict_100[mol_size_iterate] = random.sample(population=my_population, k=SAMPLE_SIZE)
        except ValueError:
            dict_100[mol_size_iterate] = my_population
            print(f'!!!\n'
                  f'Sample larger than population for molecule size: {mol_size_iterate}. '
                  f'I will save all {len(my_population)} molecules as a sample.'
                  f'Other Samples are {SAMPLE_SIZE}. So be caution!')
    # <-- ready-to-use dict with 100 or other number of molecules of every size

    # save txt files with mol idxs for every size. Folders are: 100/10, 100/11, ...
    current_path = DATA_FOLDER + '/' + str(SAMPLE_SIZE)  # i.e., folder 100 ...
    if not os.path.exists(current_path):
        os.mkdir(current_path)
    for i in MOL_SIZES_OF_INTEREST:
        file_name = current_path + '/' + str(i)
        np.savetxt(fname=file_name, X=dict_100[i], fmt='%i')

    # get smiles <-- make a separate file
    # load db
    suppl = Chem.SmilesMolSupplier(DATA_PATH)  # suppl

    # save smiles files to 100/10, 100/11, ...
    mol_100_dict = {}  # not necessary for 100 molecules. Can also be other number. But the 100 is the goal.
    smiles_100_dict = {}  #
    for atoms_number in MOL_SIZES_OF_INTEREST:
        mol_100_dict[atoms_number] = [suppl[i] for i in dict_100[atoms_number]]
        smiles_100_dict[atoms_number] = [Chem.MolToSmiles(suppl[i]) for i in dict_100[atoms_number]]
        with Chem.SmilesWriter(DATA_FOLDER + '/' + str(SAMPLE_SIZE) + '/' + str(atoms_number) + '.smi') as w:
            for m in mol_100_dict[atoms_number]:
                # print(m)
                w.write(m)

    print('I am done')


if __name__ == '__main__':
    main()
