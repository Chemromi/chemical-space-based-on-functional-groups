# Chemical Space Based on Functional Groups Toolkit

欢迎使用基于官能团的化学空间建立工具！本项目专注于向一组初始分子添加官能团来得到一小片化学空间。我们的工具包生成了大量新分子，这些分子被存储在一个CSV文件（mcsv）中，并且在另一个CSV文件（vcsv）中记录了每一次生成过程。



## Requirements

在开始之前，请确保您已经安装了以下包：

- RDKit (version > 2023)
- DGL (version > 1)
- PyTorch
- pandas



## How to Use

该工具包旨在通过命令行界面使用。为了生成新分子，您需要提供输出文件的路径，以及输入分子和官能团的列表。您还将指定生成迭代的次数。

以下是您将使用的命令行参数的简要概述：

- `--mpath`: 生成分子将被存储的CSV文件的路径。
- `--vpath`: 每一次生成过程存储的CSV文件。
- `--input`: 作为初始结构使用的输入分子列表。
- `--lines`: 要添加到初始分子的官能团列表。
- `--num`: 要执行的生成迭代次数。

### Example Command

```bash
python main.py  --mpath m1.csv --vpath v1.csv --input 'c1ccccc1' --lines 'C1=CSC=C1' 'C1=CC=NC=C1' 'CCl' 'CO' 'C1=CC=CC=C1' 'C=O' 'CN' 'CC' --num 5 &This command will generate new molecules by adding hydroxyl (OH) and amino (NH2) groups to the specified initial molecules for 100 iterations.
```

该命令将通过向指定的初始分子（苯` 'c1ccccc1'`）添加指定的官能团（噻吩、吡啶、氯甲基、甲氧基、苯基、羰基、甲胺基和甲基）进行5次迭代来生成新分子。每次迭代都将探索将每个官能团添加到不断增长的分子集中。

请将` 'c1ccccc1'` 替换为您的初始分子的SMILES，并将 `--lines` 后面的列表替换为您希望添加的官能团，也使用SMILES。根据需要的生成迭代次数调整 `--num` 参数。



## Acknowledgments

该工具的代码受到了DGL库的极大启发，特别是[dgllife.model.model_zoo.jtnn.jtnn_vae](https://github.com/awslabs/dgl-lifesci/tree/master/python/dgllife/utils/jtvae).

## License

Distributed under the Apache License 2.0. See `LICENSE` for more information.


