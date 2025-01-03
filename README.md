# Exposome
This repository will contain the analyses for the human chemical exposome
![exposome](img/exposome.png)
---
Chemical exposures exert a significant impact on human health, both on an individual and a population-wide level. To date, there is no unifying view of how diverse chemical compounds may interfere with biological processes and contribute to disease risk. In this study, we used a network-based approach to construct a comprehensive map that links 9,887 exposures through their shared impact at the genetic level. The resulting map can be used to define classes of exposures that affect the same biomolecular processes, even if they are chemically distinct. We found that exposures target specific modules within the human interactome of protein-protein interactions and that the harmfulness of an exposure is related to its interactome connectivity. A systematic comparison between the interactome modules affected by exposures and disease associated modules suggested that their proximity in the interactome can be used to predict exposure-disease relationships. To validate these predictions, we cross-referenced nation-wide disease prevalence data with reports of environmental exposures. We found that elevated levels of a particular exposure in air or water correlated with an increase in the incidence of diseases whose interactome modules overlapped with those of the respective exposures. Taken together, our study provides a blueprint for the systematic investigation of chemical exposures and their relationships, molecular impact and related disease associations, as well as population health impact.

---

### REQUIREMENTS

Create a virtual environment to install all required packages:

+ create a virtual environment
```
python3 -m venv name_of_env
```

+ activate it
```
source name_of_env/bin/activate
```

+ install requirements packages
```
python3 -m pip install -r exposome_requirements.txt
```

+ to use environment with jupyter notebooks
```
ipython kernel install --user --name=name_of_env
```

---

### First steps into the Exposome world
Please download all files from our zenodo [repository](https://zenodo.org/records/10829457) and deposit them into your directory of notebooks

```
The zip file has the following folders
├── input: these files are absolutely needed to run the analyses from scratch
├── intermediate: these files can be produced if you run all notebooks in the exact order, but computation time might be long, so we would recommend to import them
└──output: these files are the results of the analyses. You should be able to reproduce these files quite quickly by following the proper notebook order and providing both input and intermediate files.
```

---
### Exposome exploration
You can navigate directly through to the chemical projection of the Exposome, the Exposure-Exposure Network (EEN), by clicking [here](http://exposome.westeurope.cloudapp.azure.com:5000/preview), or you can connect to the dedicated [webpage](https://menchelab.com/exposomeapp). The network exploration can be performed through the link in a desktop version, or in an immersive environment with a Virtual Reality (VR) headset.


### SYSTEM REQUIREMENTS
All computation contained in the main jupyter notebook were carried out on a machine with a Apple M2 chip 8-Core CPU and 16GB of Memory.
Heavier computation, that can be found in the separate python scripts, was calculated on a cluster.
