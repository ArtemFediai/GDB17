from decimal import Decimal
from rdkit import Chem

DATA_PATH = '/home/artem/dataset_gdb17/GDB17.50000000.smi'


def main():
    print('I run main')

    suppl = Chem.SmilesMolSupplier(DATA_PATH)
    data_length = len(suppl)
    print(f'length of the data set: {Decimal(data_length):.2E}')

    for i in range(1000):
        num_mols = suppl[i].GetNumAtoms()
        print(f'num_mols: {num_mols}')

    print('I am done')

if __name__ == '__main__':
    main()
