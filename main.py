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
MOL_NUM_LIMIT = 1_000_000
# MOL_NUM_LIMIT = False  # if False, the whole data set

RDLogger.DisableLog('rdApp.*')


def main(MOL_NUM_LIMIT=None):
    print('I run main')

    suppl = Chem.SmilesMolSupplier(DATA_PATH)
    data_length = len(suppl)
    print(f'length of the data set: {Decimal(data_length):.2E}')

    if not MOL_NUM_LIMIT:
        MOL_NUM_LIMIT = data_length  # <---!
    else:
        MOL_NUM_LIMIT = MOL_NUM_LIMIT

    # randomly permuted molecule numeration
    permuted_numbers_int_64 = np.random.permutation(data_length).astype(int)
    permuted_numbers = permuted_numbers_int_64.tolist()

    numbers_of_atoms = [10, 11, 12, 13, 14, 15, 16, 17]

    idx_vs_num_atoms = {}
    for i, num_atoms in enumerate(numbers_of_atoms):
        idx_vs_num_atoms[num_atoms] = i

    atoms_count_list = [suppl[i].GetNumAtoms() for i in permuted_numbers[:MOL_NUM_LIMIT]]  # takes 1 h

    unique_atom_counts = np.unique(atoms_count_list)
    occurences_number = np.bincount(atoms_count_list)

    ax = seaborn.countplot(atoms_count_list)
    plt.show()

    # save everything
    save_dict = {
        'atom_count_list.pkl': atoms_count_list,
        'unique_atom_counts.pkl': unique_atom_counts,
        'occurences_number.pkl': occurences_number
    }

    for save_keys in save_dict:
        with open(DATA_FOLDER + '/' + save_keys, 'wb') as fid:
            pickle.dump(save_dict[save_keys], fid)


if __name__ == '__main__':
    main(MOL_NUM_LIMIT=MOL_NUM_LIMIT)
