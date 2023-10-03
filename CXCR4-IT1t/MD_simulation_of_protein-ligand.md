# How to run a Molecular Dynamics simulation of protein-ligand complex

![Include](/images/cover_figure.png)

## 1. Context of the project

In this project, we will work with a protein-ligand complex (CXCR4-IT1t). CXCR4, also known as C-X-C chemokine receptor type 4, is a G-protein coupled receptor (GPCR) that is primarily expressed on the surface of various cells, including hematopoietic stem cells, immune cells, and cancer cells. It plays a crucial role in cell migration, and hematopoiesis, among others. The natural ligand for CXCR4 is the chemokine stromal cell-derived factor-1 (SDF-1), also known as CXCL12. The CXCR4/CXCL12 complex is involved in various physiological processes and has been implicated in several pathological conditions, including cancer metastasis, HIV infection, and inflammatory diseases.<sup><a href="#ref1">[1]</a></sup>

### Importance in Drug Discovery

Given its role in various diseases, CXCR4 has emerged as an attractive target for therapeutic intervention. In cancer, the CXCR4/CXCL12 axis is involved in tumor growth, angiogenesis, metastasis, and resistance to therapy. Thus, targeting CXCR4 can potentially inhibit tumor progression and metastasis. In the context of HIV, CXCR4 acts as a co-receptor for certain strains of the virus, making it a potential target for antiretroviral therapy. Moreover, its role in inflammatory diseases has also made it a target for anti-inflammatory drugs.<sup><a href="#ref2">[2]</a></sup>

### Inhibitors and IT1t

Several small molecule inhibitors of CXCR4 have been developed for therapeutic purposes. These inhibitors can block the binding of CXCL12 to CXCR4, thereby inhibiting the downstream signaling pathways. One such inhibitor is IT1t, a small molecule that has shown potent antagonistic activity against CXCR4. IT1t has been studied for its potential in treating HIV infection and certain types of cancer. It is known to inhibit the binding of the virus to the CXCR4 co-receptor, thus preventing the entry of the virus into the host cells. In the context of cancer, IT1t can potentially inhibit tumor growth and metastasis by blocking the CXCR4/CXCL12 axis.<sup><a href="#ref3">[3]</a></sup>

### Protein-Ligand MD Simulations

When working with ligands, especially small molecules, one of the most challenging aspects is assigning an appropriate force field. This often requires using external resources or tools to ensure accurate force field parameterization. As the majority of marketed drugs are small molecules, improving the small-molecule force field parameters in ligand binding studies is crucial. Achieving improved parametrization of ligand torsion angles compared to existing force fields can significantly impact future drug design investigations.<sup><a href="#ref4">[4]</a></sup>

In this context, we will use the CGenFF (CHARMM General Force Field)<sup><a href="#ref5">[5]</a></sup>, which can generate parameters for small molecules to be used in conjunction with the CHARMM force field.

_**#NOTE:**_ CXCR4 is a membrane protein. A proper modelling of the membrane has to be handled before using it. However, the purpose of this project is to exemplify how to proceed in the design of a protein-ligand complex and, therefore, how to add a force field to a small non-protein-like molecule, which is of high importance in drug discovery. Hence, no modelling of the membrane will be performed here.

_**#NOTE 2:**_ For the sake of decreasing the calculations’ time, the timescale was chosen to be shorter than in a normal MD simulation.

## 2. System Preparation

Ensuring a well-prepared initial system is crucial for the success and accuracy of the simulation.

### <u>2.1 Prepare the Protein Topology</u>

The first step involves retrieving the CXCR4 and ligand structures. The whole structure can be found in the PDB database as 3oe6. It can be found in the folder “raw-protein”.

- #### **2.1.1 Clean the protein-ligand by removing the solvent and any other molecules.** 

You can do this by opening the file in a simple text editor and removing the water molecules (HOH in the PDB file) or by using the grep command, as exemplified here:

```
grep -v HOH protein.pdb > protein_clean.pdb
```

- #### **2.1.2 Look and recover missing residues using MODELLER.**

Before using the protein for any MD simulation, it is necessary to always inspect the .pdb file for entries labelled as MISSING RESIDUES. This means that atoms of entire parts of the protein are missing in the crystal structure. If the missing part is on the terminal region, this might not be a problem for the MD simulation. However, any internal sequences with missing atoms or amino acid residues with absent atoms will lead to errors during the creation of the topology file.

<img src="https://github.com/Stef0916/computational_chemistry/blob/46dda26d20292fbc19de558738535ab323901589/CXCR4-IT1t/images/missing%20residues.png" title="Missing Residue">

In the folder [“/protein/missing_residues”](/protein/missing_residues/) you will find all the files and scripts used to retrieve the missing residues of the protein if needed. Please, follow the instructions in MODELLER for further details.

- #### **2.1.3 Separate the protein from the ligand**
```
grep ligand protein_clean.pdb > ligand.pdb
```

The ligand and protein .pdb files can be found in the “ligand” “protein” folder as [“itd.pdb”](ligand/itd.pdb) and [“3oe6_all-residues.pdb”](protein/3oe6_all-residues.pdb), respectively.

- #### **2.1.4 Protein topology file**

I used the **CHARMM36 force field**. This can be obtained from the . In this project, the ff can be found in the “/protein/” folder under the name [“charmm36_ljpme-jul2022.ff.tgz”](/Users/stefani/Desktop/Jupyter-notebook/Molecular_Dynamics/CXCR4-IT1t/CXCR4-IT1t/ligand/charmm36_ljpme-jul2022.ff.tgz)

```
gmx pdb2gmx -f protein_clean.pdb -o protein_processed.gro -water spc
```
```
gmx: This is the main command-line interface for GROMACS.
```

- **pdb2gmx:** A tool within GROMACS that converts PDB files into GROMACS format and generates a topology.

- **-f protein_clean.pdb:** Specifies the input file.

- **-o protein_processed.gro:** Specifies the output structure file, which will be named processed.gro.

- **water spc:** This option specifies the water model to be used, which in this case is spc.

When the system asks which ff you want for your system, you’ll need to choose the CHARMM36 all-atoms force field (July 2022) from ‘/your_location/’, which usually is the first option.

This command will generate three new files:

* **processed.gro:** This is the output structure file you specified with the -o flag. It contains the coordinates of the system in GROMACS format.

* **topol.top:** This is the topology file that describes the molecular system in terms of molecules, atoms, bonds, angles, dihedrals, and other molecular parameters. It is used in subsequent steps of the simulation.

- **posre.itp:** This file contains position restraint information. Position restraints are often used during equilibration to prevent large-scale motions while allowing the system to relax.



### <u>2.2 Prepare the Ligand Topology</u>

- #### **2.2.1 Add hydrogen atoms and produce .mol2 file**

The hydrogen atoms are normally not solved during X-ray crystal structure. To do both (add H and get a .mol2 file), I used the free software AVOGADRO.

- #### **2.2.2 Correct the .mol2 file.**

Following the recommendations of [MD-tutorials](http://www.mdtutorials.com/gmx/complex/):

1. Open the .mol2 file in a text editor

2. Replace **** at the heading with the ligand’s name.

3. Fix the residue and the number of the ligand.

- #### **2.2.3 Sort bonds in ascending order**

From the previous page I mentioned, I downloaded a script (you can find it in [“/ligand/ sort_mol2_bonds.pl”](ligand/sort_mol2_bonds.pl), which helped to fix the order of the residues of the ligand.

```
perl sort_mol2_bonds.pl ligand.mol2 ligand_sorted.mol2
```


- #### **2.2.4 Generates the topol.top with CGenFF**

In the CGenFF server, after logging in:

1. Upload the ligand_sorted.mol2

2. Download the ligand_sorted.str

3. Run the py script to generate the ligand.itp and ligand.prm files needed to work with GROMACS (the script can be found in the [/ligand/ligand/cgenff_charmm2gmx_py3_nx1.py](ligand/cgenff_charmm2gmx_py3_nx1.py) folder, or in the previous cited webpage).

```
python cgenff_charmm2gmx.py LIGAND ligand_sorted.mol2 ligand.str charmm36-jul2022.ff
```


### <u>2.3 Build the complex (protein-ligand)</u>

- #### **2.3.1 Convert the ligand.pdb into .gro file**

```
gmx editconf -f ligand.pdb -o ligand.gro
```

- #### **2.3.2 In the complex.gro file, I added the coordinates of the ligand.gro at the end of complex.gro before the box vector.**

![Coordinates of the Ligand](/images/coordinates_of_the_ligand.png)


### <u>2.4 Build the complex topology file</u>

- #### **2.4.1 ligand.itp:** insert the **#include** into the topol.top after the position restrain of the protein

![Include](/images/include.png)

- #### **2.4.2 Include new geometrical parameters from the ligand.prm at the top of the topol.top, after “include forcefield parameters”**

![Include](/images/geometrical_parameters%20.png)

- #### **2.4.2 At the end of the topol.top file include the ligand and number of molecules.**

![Include](/images/ligand_and_number_of_molecules.png)


### <u>2.5 Defining the Unit Cell and Solvate the System</u>

```
gmx editconf -cp complex.gro -o complex_box.gro -bt square -d 1.0
```
```
gmx solvate -cp complex_box.gro -cs spc216.gro -o complex_solv.gro -p topol.top
```

* **solvate:** A tool within GROMACS used to immerse the molecule in a solvent.

* **editconf:** This is a GROMACS tool used to edit the structure files, such as changing the box dimensions, centring a molecule in the box, or converting file formats.

* **-cp complex.gro:** The -cp flag specifies the input structure file, which in this case is complex.gro.

* **-o complex_box.gro:** The -o flag specifies the output structure file, which will be named complex_box.gro.

* **-bt square:** The -bt flag specifies the type of box to be used. In this case, a square (cubic) box is chosen.

* **-d 1.0:** The -d flag specifies the distance (in nm) to maintain between the molecule and the edge of the box. Here, a distance of 1.0 nm is chosen, ensuring that the molecule is at least 1.0 nm away from the box boundaries.

* **-c complex_box.gro:** Specifies the input structure file, which is processed.gro.

* **-cs spc216.gro:** Specifies the solvent configuration. Here, spc216.gro is used.

* **-ocomplex_solv.gro:** Specifies the output structure file, which will be named complex_solv.gro.

* **-p topol.top:** Specifies the topology file, file to be updated.

### <u>2.6 Adding ions</u>

```
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr
```
```
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```


* **-f ions.mdp:** The -f flag specifies the input parameter file, which in this case is ions.mdp. This file contains the simulation parameters. It can be found in the [“/complex/ions.mdp”](complex/ions.mdp) folder.

* **-c complex_solv.gro:** The -c flag specifies the input structure file. This file contains the coordinates of the system.

* **-p topol.top:** The -p flag specifies the topology file to be updated.

* **-o ions.tpr:** The -o flag specifies the output run-input file, which will be named ions.tpr. This file is a compiled binary file that contains all the information needed to run the simulation.

* **genion:** A tool within GROMACS used to add ions to the system.

* **-s ions.tpr:** Specifies the input run-input file, which is ions.tpr.

* **-o solv_ions.gro:** Specifies the output structure file, which will be named solvated_ions.gro.

* **-pname NA:** Specifies the name of the positive ion, which in this case is sodium (NA).

* **-nname CL:** Specifies the name of the negative ion, which in this case is chloride (CL).

* **-neutral:** This option ensures that the system is neutralized by adding the appropriate number of ions.



## 3. Energy Minimization

Before starting the dynamics, relax the system to rectify any unfavourable contacts or configurations.

### <u>3.1 Prepare the system:</u>

```
gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
```


Execute the energy minimization:


```
gmx mdrun -v -deffnm em
```


After this step, ensure convergence by checking the energy minimization graph [/images/potential.xvg](images/potential.xvg)


![Include](/images/potential.png)


## 4. Equilibration

To stabilize the system at the desired temperature and pressure. Here are some considerations that need to be taken into account:

1. Apply restraints to the ligand:

For this, it is necessary to create an index group for the ligand containing only its non-hydrogen atoms:

```
gmx make_ndx -f ligand.gro -o index_ligand.ndx

…

> 0 & ! a H*

> q
```
```
gmx genrestr -f ligand.gro -n index_ligand.ndx -o posre_ligand.itp -fc 1000 1000 1000
```


Now, the ligand restraints need to be added to the topology file:

![Include](/images/posre_lig.png)


2. Treatment of temperature of coupling groups: 

The best is to consider the system as a single entity before the temperature equilibration step. For that, another index must be created:


```
gmx make_ndx -f em.gro – o inde.ndx

> 1 | 13

> q
```


In the nvt.md add in the tc-group = Protein_ligandname_water_and_ions



### <u>4.1 NVT Equilibration</u>

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```


* **grompp:** This is a GROMACS tool used to preprocess and compile the system for simulations. It checks the consistency of the input files, processes the topology, and compiles the system into a run-input file for molecular dynamics simulations.

* **-f nvt.mdp:** The -f flag specifies the input parameter file. This file contains the simulation parameters. It can be found in the [complex/nvt.mdp](complex/nvt.mdp) folder.

* **-c em.gro:** The -c flag specifies the input structure file, which in this case is em.gro. This file contains the coordinates of the system.

* **-r em.gro:** The -r flag specifies the reference structure for position restraints. In this case, em.gro is used as the reference structure. This is often used during equilibration phases to prevent large-scale motions of the system.

* **-p topol.top:** The -p flag specifies the topology file to be updated.

* **-o nvt.tpr:** The -o flag specifies the output run-input file, which will be named nvt.tpr. This file is a compiled binary file that contains all the information needed to run the simulation.



#### **4.1.1 Execute the NVT equilibration**

```
gmx mdrun -deffnm nvt
```


### <u>4.2 NPT Equilibration</u>

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```


#### **4.2.1 Execute the NPT equilibration**

```
gmx mdrun -deffnm npt
```


After this step, ensure convergence by checking the pressure and temperature  graph ([“/images/temperature.xvg”]() and [“/images/pressure.xvg”]()). The RMSD of the backbone is also a good practice to be checked after equilibration. The files are stored in the “/images/ folder.

![Include](/images/pressure.png) ![Include](/images/temperature.png)

After equilibration, several output files are generated. The exact files you get depend on the settings you use in your .mdp file. However, commonly generated files after equilibration include:

* **.edr:** The energy file. It contains energy data, such as potential energy, kinetic energy, temperature, pressure, and other thermodynamic properties, recorded at each saved time step.

* **.log:** The log file. It contains a detailed record of the simulation, including settings, performance data, and any warnings or errors. It's a text file and can be viewed with any text editor.

* **.gro:** This is the coordinate file that contains the final positions of all atoms in the system at the end of the equilibration. It's in GROMACS format and can be used as the starting structure for subsequent simulations.

* **.cpt:** The checkpoint file. It contains the complete state of the simulation and allows for continuation or branching of simulations. It includes information like velocities, random seeds, and more. This file is crucial if you want to continue the simulation from where it left off.

* **.xtc:** This is a compressed trajectory file that contains only the positions of the atoms (no velocities or forces). It's often used for analysis because of its smaller size compared to the .trr file.

* **.tpr:** It compiles the system’s initial atomic coordinates, simulation parameters (from the .mdp file), force field parameters, box dimensions and type. All of which is needed to run a simulation.



## 5. Production



The production phase generates trajectories and data for subsequent analysis. Now, the restraints are released, and the system is ready to run the simulation.



### 5.1 <u>Prepare the system:</u>

```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_2.tpr
```


* **grompp:** This is a GROMACS tool used to preprocess and compile the system for simulations. It checks the consistency of the input files, processes the topology, and compiles the system into a run-input file for molecular dynamics simulations.

* **-f md.mdp:** The -f flag specifies the input parameter file, which in this case is md.mdp. This file contains the simulation parameters. It can found in [complex/md.mdp](complex/md.mdp)

* **-c npt.gro:** The -c flag specifies the input structure file, which in this case is npt.gro. This file contains the coordinates of the system.

* **-t npt.cpt:** The -t flag specifies the input checkpoint file, which in this case is npt.cpt. The checkpoint file contains the state of the previous simulation, allowing for continuation or branching of simulations. It includes information like velocities, random seeds, and more.

* **-p topol.top:** The -p flag specifies the topology file to be updated.

* **-n index.ndx:** The -n flag specifies the index file,. The index file contains groups of atoms that can be used to define specific parts of the system for analysis or other operations. In this case, it is the index of the protein and the ligand.

* **-o md_0_2.tpr:** The -o flag specifies the output run-input file. This file is a compiled binary file that contains all the information needed to run the simulation.



### 5.2 <u>Execute the production MD run:</u>

```
gmx mdrun -deffnm md_0_02
```


## 6. Analysis

Bring protein in the center of the box

```
gmx trjconv -s md_0_02.tpr -f md_0_02.xtc -o md_0_02_center.xtc -center -pbc mol -ur compact
```


After this, depending on the purpose of the MD simulation, a huge set of analyses can be performed. In this particular case, I was interested for example, about the interaction between the ligand and the protein residues. 

![Include](/images/it1t-protein_1_distance.png)
![Include](/images/it1t-protein_2_distance.png)
![Include](/images/it1t-protein_3_distance.png)
![Include](/images/it1t-protein_4_distance.png)


## 7. References

<ol>
<li id="ref1"><a href="#" target="_blank">Teicher, B. A., & Fricker, S. P. (2010). CXCL12 (SDF-1)/CXCR4 pathway in cancer. Clinical Cancer Research, 16(11), 2927-2931.</a></li>
<li id="ref2"><a href="#" target="_blank">Balkwill, F. (2004). The significance of cancer cell expression of the chemokine receptor CXCR4. Seminars in Cancer Biology, 14(3), 171-179.</a></li>
<li id="ref3"><a href="#" target="_blank">De Clercq, E. (2003). The bicyclam AMD3100 story. Nature Reviews Drug Discovery, 2(7), 581-587.</a></li>
<li id="ref4"><a href="#" target="_blank">Raniolo, S. & Limongelli, V. (2021). Improving Small-Molecule Force Field Parameters in Ligand Binding Studies. Frontiers in Molecular Biosciences, 8.</a></li>
<li id="ref5"><a href="#" target="_blank">Vanommeslaeghe, K., Hatcher, E., Acharya, C., Kundu, S., Zhong, S., Shim, J., Darian, E., Guvench, O., Lopes, P. E. M., Vorobyov, I., & MacKerell, A. D. (2009). CHARMM general force field: A force field for drug‐like molecules compatible with the CHARMM all‐atom additive biological force fields. Journal of Computational Chemistry, 31.</a></li>
</ol>