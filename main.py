import os.path
from decimal import Decimal
import numpy as np
import pickle
# plot
import matplotlib.pyplot as plt
import seaborn
# chemistry
from rdkit import Chem
from rdkit import RDLogger

DATA_FOLDER = '/home/artem/dataset_gdb17'
DATA_FILE_NAME = 'GDB17.50000000.smi'
MOL_NUM_LIMIT = 1_000_000  # how many molecules from the data set
# MOL_NUM_LIMIT = False  # if False, the whole data set

# warning
RDLogger.DisableLog('rdApp.*')  # disable warning from rdkit


def main(MOL_NUM_LIMIT=None):
    print('I run main')
    DATA_PATH = DATA_FOLDER + '/' + DATA_FILE_NAME
    suppl = Chem.SmilesMolSupplier(DATA_PATH)
    data_length = len(suppl)
    print(f'length of the whole data set: {Decimal(data_length):.2E}')

    if not MOL_NUM_LIMIT:
        MOL_NUM_LIMIT = data_length  # <---!
    else:
        MOL_NUM_LIMIT = MOL_NUM_LIMIT

    # randomly permuted molecule numeration
    permuted_numbers_int_64 = np.random.permutation(data_length).astype(int)  # 64bit does not work for rdkit
    permuted_numbers = permuted_numbers_int_64.tolist()

    # numbers_of_atoms = [10, 11, 12, 13, 14, 15, 16, 17]

    # idx_vs_num_atoms = {}
    # for i, num_atoms in enumerate(numbers_of_atoms):
    #     idx_vs_num_atoms[num_atoms] = i

    atoms_count_list = [suppl[i].GetNumAtoms() for i in permuted_numbers[:MOL_NUM_LIMIT]]  # takes 1 h

    unique_atom_counts = np.unique(atoms_count_list)
    occurences_number = np.bincount(atoms_count_list)

    ax = seaborn.countplot(atoms_count_list)
    fig = ax.get_figure()

    # save plot num occurences of X-atoms molecule
    current_path = DATA_FOLDER + '/' + 'atoms_count_list.png'
    if not os.path.exists(DATA_FOLDER + '/' + 'atoms_count_list.png'):
        print(f'Saving atoms_count_list into {current_path}...')
        fig.savefig(current_path)
        print('...Saved successfully')
    else:
        print(f'the path {current_path} exists. Will not save the Figure')

    # save data
    save_dict = {
        'atom_count_list.pkl': atoms_count_list,
        'unique_atom_counts.pkl': unique_atom_counts,
        'occurences_number.pkl': occurences_number,
        'permuted_numbers.pkl': permuted_numbers,
    }

    for save_key in save_dict:
        current_path = DATA_FOLDER + '/' + save_key
        if not os.path.exists(current_path):
            print(f'Saving {save_key} in {current_path}...')
            with open(current_path, 'wb') as fid:
                pickle.dump(save_dict[save_key], fid)
            print(f'...Saved successfully')
        else:
            print(f'the path {current_path} exists. Will not save anything')


if __name__ == '__main__':
    main(MOL_NUM_LIMIT=MOL_NUM_LIMIT)
