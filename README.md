# Network model to study cell division orientation in epithelia monolayers

To gain a conceptual understanding of how tension in mitotic cells evolves relative to the tension 
in interphase cells and in the tissue, we devised a simple computational model of the monolayer as 
a purely elastic 2D material in which cells have a preset active tension and rigidity. By varying tension, 
rigidity and position of the boundaries, we aim to reproduce a range of experimental conditions and 
characterize the distribution of stress in the viscinity of the mitotic cell. 
The simulated scenarios are:
1. monolayer in its initial configuration clamped at both ends is subjected to a tensile stress due to cell contractility
2. the monoalyer with internal cell contractility is compressed
3. the internal contractility is increased while the monolayer is under a compressive state
4. the internal contractility of a clamped monolayer is reduced
5. the monolayer with reduced contractility is stretched

### How to run the code

- Install Julia (latest tested version 1.8.2)
- Set the parameters in the "Main" section of the main.jl file (e.g. n: number of steps for the simulation - if one of the parameter varies; param1: EA, where E is the Young's modulus and A the cross section area; param2: displacement of the two ends, param3: prestress) 
- Change filename and figname with the directory where the images and tensors .txt files will be saved
- Run the main.jl file in REPL as include("main.jl")

## Scientific publication:

 "Tension at intercellular junctions is necessary for accurate orientation of cell division in the epithelium plane"
Ana Lisica, Jonathan Fouchard, Manasi Kelkar, Tom P. J. Wyatt, Julia Duque, Anne-Betty Ndiaye, Alessandra Bonfanti, Buzz Baum, Alexandre J. Kabla, Guillaume T. Charras

https://www.biorxiv.org/content/10.1101/2022.01.30.478396v1.full
