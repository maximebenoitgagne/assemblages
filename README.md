Welcome! This repository contains the model and the scripts used to generate the figures in
Benoît-Gagné et al.,
"Interannual controls of the phytoplankton assemblages in Baffin Bay".

Directory structure:

- The model itself is contained in the following directories: doc, eesupp, jobs, lsopt, model, optim, pkg, tools, utils and verification.
- The configuration of the model is in the directory gud_1d_35+16.
- The jupyter notebook assemblages.ipynb generates the figures in the directory
figures_assemblages_progress.
- The jupyter notebook assemblages_supmat.ipynb generates the figures in the
directory figures_assemblages_sup_progress.
- Data for the generation of the figures (Data S1 to S5). It includes complete model
outputs for the the reference experiment (EXP-0) and for 10 of the 121 sensitivity experiments. Only statistics are provided for the other 111 sensitivity experiments. The complete model outputs for all experiments are available on the Federated Research Data Repository (FRDR, see paper).

Model:

- The model is a one-dimensional configuration of the
biogeochemical/ecosystem model of Dutkiewicz et al. (2021) in
*Glob. change biol.* [https://doi.org/10.1111/gcb.15493](https://doi.org/10.1111/gcb.15493)
The tracers are mixed by the MIT general circulation model
(MITgcm, Marshall et al., 1997 in *JGR*).
The paper Benoît-Gagné et al. (2024) in *Elem. Sci. Anth.* [https://doi.org/10.1525/elementa.2024.00008](https://doi.org/10.1525/elementa.2024.00008) describes some modifications relative to
Dutkiewicz et al. (2021).

Datasets:

- data/DataS1_winter_nutrients_in_situ_qik:
Forcing fields of nitrate and silicate concentrations from January 1 to May 15 for the reference simulation (EXP-0). They are also the in situ nitrate and silicate concentrations at the Qikiqtarjuaq sea ice camps between mid-April and the end of May in 2015 and 2016 (67.4797°N, -63.7895°E). It contains a subset of the files available in the dataset Massicotte et al. (2019). [https://doi.org/10.17882/59892](https://doi.org/10.17882/59892). The paper related to this dataset is Massicotte et al. (2020). [https://doi.org/10.5194/essd-12-151-2020](https://doi.org/10.5194/essd-12-151-2020). Details are available in MetadataS1.pdf.
- data/DataS2_winter_nutrients_simulated:
Forcing fields of nitrate and silicate concentrations from January 1 to May 15 for the 121 sensitivity experiments. They don't correspond to observations. They were constructed by multiplying the observed vertical profiles by a factor to get a nutrient concentration of xx uM at 100 m where xx is in the name of the file. See [https://github.com/maximebenoitgagne/gud_groups/blob/gud/gud_1d_35%2B16/input_noradtrans/input/data_rbcs.ipynb](https://github.com/maximebenoitgagne/gud_groups/blob/gud/gud_1d_35%2B16/input_noradtrans/input/data_rbcs.ipynb) for more details. Details are available in MetadataS2.pdf.
- data/Data3_outputs_short:
Annual average of biomass and annually integrated relative contribution of each phytoplankton group for the 121 sensitivity simulations. Details are available in MetadataS3.pdf.
- data/DataS4_outputs_long:
Simulation outputs for the reference simulation (EXP-0) and for the 10 sensitivity simulations discussed in more detail in the paper (EXP-S02N02, EXP-S02N08, EXP-S02N16, EXP-S08N02, EXP-S08N08, EXP-S08N16, EXP-S14N02, EXP-S14N08, EXP-S14N16 and EXP-S14N20). Details are available in MetadataS4.pdf.
- data/DataS5_winter_nutrients_in_situ_11exp:
Winter nitrate and silicate concentrations for 11 expeditions from literature (Ardyna et al., 2020 in *Elem. Sci. Anth.* [https:doi.org/10.1525/elementa.430](https:doi.org/10.1525/elementa.430)). The csv file was provided by Mathieu Ardyna. Details are available in MetadataS5.pdf.

Notes:

Some data files are larger than 100 GB.
Hence, if Git LFS is not installed, the files larger than 100 GB will be replaced with placeholders after cloning the project.

The exact procedure I used to deploy the code on a supercomputer of the Digital Research Alliance of Canada with the SLURM workload manager is

```
module load git-lfs
git_lfs clone git@github.com:maximebenoitgagne/gud_groups.git
```

The procedure to run the model is in the README of the directory gud_1d_35+16.

Let me know if you have any requests or comments.
You can contact me via
[ResearchGate](https://www.researchgate.net/profile/Maxime-Benoit-Gagne).
