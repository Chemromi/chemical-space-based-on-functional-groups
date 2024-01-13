# Chemical Space Based on Functional Groups Toolkit

Welcome to the Chemical Space Based on Functional Groups Toolkit! This project focuses on expanding the chemical space of molecules by systematically adding functional groups to a set of initial molecules. Our toolkit generates a plethora of novel molecules, which are stored in a CSV file (mcsv), and each generation process is recorded in another CSV file (vcsv).

[中文](https://github.com/Chemromi/chemical-space-based-on-functional-groups/blob/main/README_zh.md)

## Requirements

Before you begin, make sure you have the following packages installed:

- RDKit (version > 2023)
- DGL (version > 1)
- PyTorch
- pandas

Please follow the official documentation to install these dependencies if you haven't already.

## How to Use

This toolkit is designed to be used with a command-line interface. To generate new molecules, you need to provide the paths to output files, as well as a list of input molecules and functional groups. You will also specify the number of generation iterations.

Here's a brief overview of the command-line arguments you'll be using:

- `--mpath`: Path to the CSV file where the generated molecules will be stored.
- `--vpath`: Path to the CSV file where the vectors of molecules differing by one functional group will be stored.
- `--input`: List of input molecules to be used as initial structures.
- `--lines`: List of functional groups to be added to the initial molecules.
- `--num`: The number of generation iterations to perform.

### Example Command

```bash
python main.py  --mpath m1.csv --vpath v1.csv --input 'c1ccccc1' --lines 'C1=CSC=C1' 'C1=CC=NC=C1' 'CCl' 'CO' 'C1=CC=CC=C1' 'C=O' 'CN' 'CC' --num 5 &This command will generate new molecules by adding hydroxyl (OH) and amino (NH2) groups to the specified initial molecules for 100 iterations.
```

This command will generate new molecules by adding the specified functional groups (thiophene, pyridine, chloromethyl, methoxy, phenyl, carbonyl, methylamine, and methyl groups) to benzene (`'c1ccccc1'`) for 5 iterations. Each iteration will explore the addition of each functional group to the growing set of molecules.

Please replace `'c1ccccc1'` with your initial molecule in SMILES notation, and the list after `--lines` with the functional groups you wish to add, also in SMILES notation. Adjust the `--num` parameter to the desired number of generation iterations.

![image](https://github.com/Chemromi/chemical-space-based-on-functional-groups/blob/main/image.JPG)

## Acknowledgments

The code of this tool is heavily inspired by the DGL library, especially [dgllife.model.model_zoo.jtnn.jtnn_vae](https://github.com/awslabs/dgl-lifesci/tree/master/python/dgllife/utils/jtvae).

## License

Distributed under the Apache License 2.0. See `LICENSE` for more information.


