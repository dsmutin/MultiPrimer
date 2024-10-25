# MultiPrimer
Tool for obtaining best primers based on several genomes using evolutionary algorithm

## Algorithm

- Select file which further will be used for primer generation
- Select files for checking universality of the primer
- Select files for checking specificity of the primer
- Select layouts
- Run (wrapped):
  - primer3 generation
  - blastn check
  - parsing
  - mutation in primers

## How To
See <a href = "exec.sh">example</a>
- data selection
- blast db build
- pipeline run

For blastn db build, you can use <a href = "scripts/prep_db.sh">prep_db.sh</a>

## Usage
Run `python pipeline.py`

```
Generation of primers based on fasta-files and blastn databases. To use it, select one
reference file to generate the initial primer set; blastn base to check primer universality
and cut off multimapping; blastn bases to remove non-specific primers Requires primer3 and
blastn pre-installed

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file for generation. Primers are generated for different
                        contigs separatly. Only gene-coding regions recommended (.fna)
  -tb TRUE_BASE, --true_base TRUE_BASE
                        Input blastn database path for primer adjusting
  -fb [FALSE_BASE ...], --false_base [FALSE_BASE ...]
                        Input blastn database path for non-specific testing. Wildcards are
                        not accepted
  -c CONTIG_TABLE, --contig_table CONTIG_TABLE
                        .tsv table with blast db information
  -o OUTPUT, --output OUTPUT
                        Output path
  -t THREADS, --threads THREADS
                        number of threads
  -ot OUTPUT_TMP, --output_tmp OUTPUT_TMP
                        Output .tmp dicrectory path for calculations and data processing.
                        .tmp in output directory as default
  -N ITERATIONS, --iterations ITERATIONS
                        Maximum iterations of evolutionary algorithm. 100 by default
  -T TOP, --top TOP     Top primers to mutate and use in next generation
  -M MUTATION_RATE, --mutation_rate MUTATION_RATE
                        Mutation probability per position of primer
  -S SET_SIZE, --set_size SET_SIZE
                        Size of mutated primers per primer
  -A APPEND, --append APPEND
                        Append best primers to array in evolutionary algoritm
  --primer3 PRIMER3     primer3_core path or command to exec. 'primer3' as default
  --blastn BLASTN       blastn path or command to exec. 'blastn' as default
  --add_set [ADD_SET ...]
                        file to set of primers to append to initial primer3 generation. empty
                        by default
  --PRIMER_PICK_PRIMER PRIMER_PICK_PRIMER
                        primer3 template option. Number of primers to pick
  --PRIMER_NUM_RETURN PRIMER_NUM_RETURN
                        primer3 template option. initial set size per gene
  --PRIMER_OPT_SIZE PRIMER_OPT_SIZE
                        primer3 template option
  --PRIMER_MIN_SIZE PRIMER_MIN_SIZE
                        primer3 template option
  --PRIMER_MAX_SIZE PRIMER_MAX_SIZE
                        primer3 template option
  --PRIMER_PRODUCT_SIZE_RANGE PRIMER_PRODUCT_SIZE_RANGE
                        primer3 template option. 2 values sepatated by '-'
  --word_size WORD_SIZE
                        blastn template option
  --reward REWARD       blastn template option
  --penalty PENALTY     blastn template option
  --gapopen GAPOPEN     blastn template option
  --gapextend GAPEXTEND
                        blastn template option
  --evalue EVALUE       blastn template option
  --max_mismatch MAX_MISMATCH
                        primer_check template option. maximum avialable mismatch
  --multimap_max MULTIMAP_MAX
                        primer_check template option. maximum multimapped hits
  --negative_max NEGATIVE_MAX
                        primer_check template option. maximum negative hits
  --min_ident MIN_IDENT
                        primer_check template option. minimal identity, percent
```

