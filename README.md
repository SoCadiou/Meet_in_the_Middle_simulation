# Meet in the Middle simulation
Simulation script from Cadiou et al., 2021 (link to add after publication).
Instructions on how to use Rscripts are available on the Rscripts themselves.
Seeds need to be specified in order to make results reproducible.

# General information
Scripts were developed to perform a three-layer simulation. Using real data from two biological layers (eg exposome and methylome) in a specific population sample of size p, two layers with realistic correlation are generated. 

In a first simulation script, linear relations between the two layers are added to simulate causal relationships and to generate a Gaussian outcome variable linearly dependant of the two layers (simulation_1_causal_effect). Various parameters are available to vary the effect range and the number of predictors from each layer and to choose which causal links to add.

In a second script (simulation_2_reverse_causal_effect), the outcome layer is generated from a third real dataset of dimension 1xp and linear relations in order to simulate a reverse causal link from the outcome on the two other layers can be generated.

8 methods aiming at selecting true predictors from the first (eg exposome) biological layer are implemented and their performance assessed in term of sensitivity and false detection proportion. 4 of them use the intermediary (eg methylome) layer.

In (Cadiou et al., 2021), real data from the Helix project [(Vrihjeid et al., 2014)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4048258/) were used. These data are not publicly available and thus not provided here but can be obtained under request to the [Helix consortium](https://www.projecthelix.eu/), after signature of a Data Transfer Agreement. 

For more details, please refer to (Cadiou et al., 2021).


