# DNA Sequence K-mer-Based Search for Microbiome Identification

We implemented a K-mer search for DNA sequences and analyzed its effectiveness comparing to a full alignment search in identifying microbiome DNA strands in Pittsburgh rivers.


## DNA Preparation
Water samples were collected at six different locations (including pure water as control). DNAs were extracted and sequenced beforehand in the wet-lab. 

The sequence reads are available in `data/Fall2018CleanReads.fa`. DNA reads quality control data and reports are available in `dna_quality_control/`. 

We also used `bacterial_16s_genes.fa`, which contains bacterial 16s genes as the database that we search against.

## Usage
Simply run `microbiome_test.py` and bacterial phylum plots will be generated in the `fraction_plots/' directory with our pre-tuned parameters for k-mer search ([see details here](https://en.wikipedia.org/wiki/K-mer)).

![Phylum Fractions at Neville Island Sample Point]()

An analysis report on parameter search has been provided in `kmer_size_acc_reports/` with accuracy data included. 

![Optimal: K=8, thresh=.6 -> Acc 0.973]()

We have run a full alignment (consult `alignment.py` for its implementation details) on the DNA reads against the s16 database on Google Cloud and the result is saved in `data/*.txt`. Runtime generated calculations have been saved in `cache/*.pickle`. You may perform further analysis on these.


## Contributors
Anupam Pokharel and Dr. Kangas

## License & Rules
[MIT](https://choosealicense.com/licenses/mit/)
[CMU Academic Policies](https://www.cmu.edu/policies/student-and-student-life/academic-integrity.html)

