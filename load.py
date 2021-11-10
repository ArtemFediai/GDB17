from decimal import Decimal
from rdkit import Chem
from rdkit import RDLogger
import numpy as np
import matplotlib.pyplot as plt
import seaborn
import pickle

DATA_FOLDER = '/home/artem/dataset_gdb17'
DATA_FILE_NAME = 'GDB17.50000000.smi'
DATA_PATH = DATA_FOLDER + '/' + DATA_FILE_NAME
# MOL_NUM_LIMIT = False  # if False, the whole data set
RDLogger.DisableLog('rdApp.*')


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
    occurences_number = save_dict['occurences_number.pkl']
    permuted_numbers = save_dict['permuted_numbers.pkl']

    ax = seaborn.countplot(atoms_count_list)
    plt.show()

    print(f'Number of atoms per molecules in the dataset: {unique_atom_counts}')

if __name__ == '__main__':
    main()
