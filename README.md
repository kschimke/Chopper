# Chopper 
## Assembling Plasmid Consensus Reads

Chopper takes consensus reads (C3POa / TideHunter / similar tool output) of ONT plasmid sequencing data, determines plasmid size(s), and generates assemblies for individual plasmids. 

--------------------------------------------------------------------------------

## Dependencies

- [Python 3](https://www.python.org/downloads/)
- [NumPy](https://pypi.org/project/numpy/)
- [minimap2](https://github.com/lh3/minimap2)
- [mappy](https://pypi.org/project/mappy/)
- [samtools](https://www.htslib.org/download/)
- [racon](https://github.com/isovic/racon)
- [medaka](https://github.com/nanoporetech/medaka/blob/master/docs/installation.rst)

--------------------------------------------------------------------------------

## Usage

```bash
python3 Chopper.py mode 
                -c consensus.fasta
                -s subreads.fastq
                -o /path/to/output/directory/ 
                -b bin size
                -p plasmid size
                -t assembly threshold
                -su subsample 
```

Arguments:
```
mode Comma separated list of modes to run. Defaults to detect,assemble
     detect:
     assemble:
     plot:

-c   consenus reads in FASTA format

-s   subreads in FASTQ format, optional if not using C3POa output

-o   output path, defaults to current working directory 

-b   bin size when sorting reads by length and used to define range around user defined plasmid size

-p   precise plasmid size, or comma separated list of sizes  

-r   range of bins to print in plot mode

-t   minimum percentage of total consensus reads required to consider a bin a peak in detect mode

-su  maximum number of reads used to make an assembly

-V   verbose, prints tool output to console otherwise creates log files

-v   print the Chopper version and exit

-h   print this help text and exit
```
