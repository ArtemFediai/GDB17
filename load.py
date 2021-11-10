# load and get 100 mols of every atoms count
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
SAMPLE_SIZE = 100

def main():
    print('I run main')

    atoms_count_list = []
    unique_atom_counts = []
    occurences_number = []
    permuted_numbers = []

    save_dict = {
        'atom_count_list.pkl': atoms_count_list,
        'unique_atom_counts.pkl': unique_atom_counts,
        'occurences_number.pkl': occurences_number,
        'permuted_numbers.pkl': permuted_numbers
    }

    for save_keys in save_dict:
        with open(DATA_FOLDER + '/' + save_keys, 'rb') as fid:
            save_dict[save_keys] = pickle.load(fid)

    atoms_count_list = save_dict['atom_count_list.pkl'] #
    unique_atom_counts = save_dict['unique_atom_counts.pkl']  # number of atoms
    occurences_number = save_dict['occurences_number.pkl']  # 0 atoms, 1 atom, 2 atoms, ...
    permuted_numbers = save_dict['permuted_numbers.pkl']

    ax = seaborn.countplot(atoms_count_list)
    plt.show()

    # print(f'Number of atoms per molecules in the dataset: {unique_atom_counts}')
    # print(f'Count of molecules with respecting atom number: {occurences_number}')

    dict_atoms_num_mol_count = dict(zip(np.arange(len(occurences_number)), occurences_number))
    print(dict_atoms_num_mol_count)

    dict_num_atoms_mol_idx = {}
    for mol_size_iterate in MOL_SIZES_OF_INTEREST:
        print(f'mol size of interest = {mol_size_iterate}')
        dict_num_atoms_mol_idx[mol_size_iterate] = [mol_idx for mol_idx, mol_size in enumerate(atoms_count_list) if mol_size == mol_size_iterate]

    current_path = DATA_FOLDER + '/' + 'dict_num_atoms_mol_idx.pkl'
    with open(current_path, 'wb') as fid:
        pickle.dump(dict_num_atoms_mol_idx, fid)

    # dict that contains mol idx of random molecules with X number of atoms. X = 10...18 (why 18 is included, no idea)
    dict_100 = {}
    for mol_size_iterate in MOL_SIZES_OF_INTEREST:
        dict_100[mol_size_iterate] = random.sample(population=dict_num_atoms_mol_idx[mol_size_iterate], k=SAMPLE_SIZE)

    current_path = DATA_FOLDER + '/' + str(SAMPLE_SIZE)
    if not os.path.exists(current_path):
        os.mkdir(current_path)
        for i in MOL_SIZES_OF_INTEREST:
            file_name = current_path + '/' + str(i)
            np.savetxt(fname=file_name, X=dict_100[i], fmt=)

    print('I am done')
if __name__ == '__main__':
    main()
