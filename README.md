
Development of Multiscale Ultra-Coarse-Grained Models for the SARS-CoV-2 Virion from Cryo-Electron Microscopy Data
=======

Authors
-------

Fengyu Li, Yuwei Zhang, Fei Xia, and Xin Xu

Introduction
------------

These files are the Electronic Supplementary Information of article entitled **Development of Multiscale Ultra-Coarse-Grained Models for the SARS-CoV-2 Virion from Cryo-Electron Microscopy Data**. One can run the molecular dynamics simulations of SARS-CoV-2 virions using the models presented here: the intact ultra-coarse-grained (UCG) model, the UCG model without ribonucleoproteins (RNPs), and the all-atom(AA)/UCG model. The AA model is extracted from `ACS Cent. Sci., 2020, 6, 1722-1734`.

Containing Files
----------------

```
.:
AA_UCG  fix  readme.md  UCG_intact  UCG_noRNP

./AA_UCG:
AA_UCG.data.gz  AA_UCG.mp4  file1_xaa.txt.gz  file2.txt
AA_UCG.in       AA_UCG.pdb  file1_xab.txt.gz

./fix:
fix_reduce.cpp  fix_reduce.h

./UCG_intact:
file1_xaa.txt.gz  file2.txt           UCG_intact.in   UCG_intact.pdb
file1_xab.txt.gz  UCG_intact.data.gz  UCG_intact.mp4

./UCG_noRNP:
file1_xaa.txt.gz  file2.txt          UCG_noRNP.in   UCG_noRNP.pdb
file1_xab.txt.gz  UCG_noRNP.data.gz  UCG_noRNP.mp4
```

File Descriptions
-----------------

|     File Types       |                Descriptions               |                  Notes                 |
| :------------------: | :---------------------------------------- | :------------------------------------- |
|      `.data.gz`      | LAMMPS data files                         | Use `gzip` to decompress `.gz` files   |
|        `.in`         | LAMMPS input files (as examples)          |                                        |
|        `.mp4`        | Movies of MD trajectories                 |                                        |
|        `.pdb`        | Structure files of models                 |                                        |
|   `.txt`, `.txt.gz`  | Parameter files of harmonic bonds         | Use `gzip` to decompress `.gz` files   |
|     `.cpp`, `.h`     | Extra `fix` style files needed for LAMMPS | Compile them before running simulation |

Simulation setup
----------------

+ Install necessary packages in LAMMPS, like `MOLECULE`, `MISC`. Then compile `fix_reduce.cpp` and `fix_reduce.h` files in LAMMPS `src/` directory:

	```sh
	make serial
	```

+ In each subdirectories, decompress the `.gz` files, and merge two `.txt` files into `file1.txt`:

	```sh
	gzip -d *.gz
	cat file1_xaa.txt file1_xab.txt > file1.txt
	```

+ Run the simulation with the sample input files. For example:

	```sh
	lmp_serial -i AA_UCG.in
	```

+ Each simulation should generate one `.dcd` trajectory file and one `log.lammps` output file.

Contact
-------

fxia@chem.ecnu.edu.cn

xxchem@fudan.edu.cn

Mar. 2023
