# Non Readable Range Finder
This repo houses a collection of scripts to find duplicate and therefore non-readable ranges in a transcriptome. The `data/gene_correction_map.csv` file contains the number of CDS indices for 

## Dependencies
- Python for taking a Fasta file and converting it into a JSON map of genes to sequence strings. RiboPy is also needed.
- Go for taking the JSON map and creating a map of duplicates.

## Input files
- A ribo file from the target organism. This is just used to get a list of genes.
- An input Fasta file for the target organism.

There are some defaults for this already set that should work on Mozart.

## Finding duplicate sequences
First, run:
```
python fasta_to_json.py
```

This will dump a `sequence_dict.json` file in the `data` directory. You can customize the ribo file and input Fasta file with optional command line arguments, but by default it will run over human genes.

Next, run the Go file. Here is the syntax:
```
go run dupe_seq_finder.go -min <min read length> -max <max read length>
```

This script is computationally intensive in the first couple minutes of execution and also memory intensive. Each execution of the script takes around 4-6 minutes. In my testing on Mozart it took up around 15GB of memory per read length analyzed, so run it in chunks as to not run out of memory. The attatched bash script can automate this. Run it with:
```
bash runner.bash
```
However, depending on what else is running on Mozart, you might be able to run it in bigger chunks. I recommend processing read lengths [21, 40].

This script will dump JSON files into the `data/dupe_seq` directory for each read length. The files are essentially giant dictionaries, where the keys represent a sequence substring and the value represents the number of times the substring appears in the transcriptome.

## Results at the gene level
Now that we have an idea of which nucleotide sequences are duplicates, we can go through all the genes again and find which indices for each gene will be non readable. We can do that by running the following script.

```
go run gene_dupe_search.go -min 21 -max 40 -filter_cds_range
```

This script will output a list of duplicate indices for each gene for each read length into `data/dupe_idx_maps.json`. It takes around 30 seconds to execute on Mozart.

Note: Make sure for whatever min and max read length you input into this script, the relevant files for those read lengths created by the previous script exist. You can check in `data/dupe_seq` to verify.

## Analyzing the output
`non_readable_range_analysis.ipynb` creates a CSV that summarizes the data (the number of non readable indices for each gene and read length) and contains some tools to analyze the data, including generation of some charts. To run this you will need the general dependencies to run a Jupyter notebook (ipykernel, matplotlib, etc.) as well as RiboPy.

