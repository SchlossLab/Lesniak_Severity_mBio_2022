## The gut bacterial community potentiates *Clostridioides difficile* infection severity.

### Abstract

The severity of *Clostridioides difficile* infections (CDI) has increased over the last few decades. Patient age, white blood cell count, creatinine levels as well as *C. difficile* ribotype and toxin genes have been associated with disease severity. However, it is unclear whether specific members of the gut microbiota associate with variation in disease severity. The gut microbiota is known to interact with *C. difficile* during infection. Perturbations to the gut microbiota are necessary for *C. difficile* to colonize the gut. The gut microbiota can inhibit *C. difficile* colonization through bile acid metabolism, nutrient consumption and bacteriocin production. Here we sought to demonstrate that members of the gut bacterial communities can also contribute to disease severity. We derived diverse gut communities by colonizing germ-free mice with different human fecal communities. The mice were then infected with a single *C. difficile* ribotype 027 clinical isolate which resulted in moribundity and histopathologic differences. The variation in severity was associated with the human fecal community that the mice received. Generally, bacterial populations with pathogenic potential, such as *Enterococcus*, *Helicobacter*, and *Klebsiella*, were associated with more severe outcomes. Bacterial groups associated with fiber degradation and bile acid metabolism, such as *Anaerotignum*, *Blautia*, *Lactonifactor*, and *Monoglobus*, were associated with less severe outcomes. These data indicate that, in addition to the host and *C. difficile* subtype, populations of gut bacteria can influence CDI disease severity.

### Importance

*Clostridioides difficile* colonization can be asymptomatic or develop into an infection, ranging in severity from mild diarrhea to toxic megacolon, sepsis, and death. Models that predict severity and guide treatment decisions are based on clinical factors and *C. difficile* characteristics. Although the gut microbiome plays a role in protecting against CDI, its effect on CDI disease severity is unclear and has not been incorporated into disease severity models. We demonstrated that variation in the microbiome of mice colonized with human feces yielded a range of disease outcomes. These results revealed groups of bacteria associated with both severe and mild *C. difficile* infection outcomes. Gut bacterial community data from patients with CDI could improve our ability to identify patients at risk of developing more severe disease and improve interventions which target *C. difficile* and the gut bacteria to reduce host damage.






### Overview

	project
	|- README         # the top level description of content (this doc)
	|- CONTRIBUTING   # instructions for how to contribute to your project
	|- LICENSE        # the license for this project
	|
	|- submission/	  # files for manuscript submission
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created
	|
	|- code/          # any programmatic code
	| |- ml/          # R scripts to run machine learning model and process data
	| |- mothur/      # dir for local copy of mothur
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make should be located in the user's PATH
* mothur (v1.XX.0) should be located in the user's PATH
* R (v. 3.X.X) should be located in the user's PATH
* etc.


#### Running analysis

```
git clone https://github.com/SchlossLab/Lesniak_Severity_XXXX_2022.git
make write.paper
```

