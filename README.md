# Code necessary to reproduce results from the BASSGWAS manuscript
_Note: to access the actively maintained version of BASSGWAS for your own experiments visit https://github.com/dhelekal/BASSGWAS_

## Setup
1. Use `Pkg.develop` to install the `julia` packages `BASSGWAS/BASSGWAS` and `BASSGWAS_old/BASSGWA_old`
2. Download the input data files from the data repository and move them into the `data` directory

## Layout
- `BASSGWAS/` Package version with varying effect scale. Used for benchmarking.
- `BASSGWAS_old/` Package version with fixed effect scale. Used for some benchmarking and lab application.
- `data/` Data directory for benchmarks. Data available in zenodo data repository.
- `run_benchmarks/` Benchmarking and post-processing launch scripts
    - `make_*.sh` Scripts to generate initial seed batches for the respective benchmarks. No needto run as outputs are included in `out`
    - `run_cip[4|8|16].sh` Scripts to launch CIP benchmarks with different batch sizes usingvarying scale
    - `run_tet[2.5|10].sh` Scripts to launch TET benchmarks at 2.5% and 10% frequency withvarying scale
    - `run_azi_mic2.sh` Script to launch AZI benchmarks with varying scale
    - `run_[azi|cip|tet]_fixedsig.sh` Scripts to launch adaptive design validation for fixedscale model
    - `annotate_all_patts.sh` Script to annotate all unitigs. Requires assemblies. No need torun, outputs are already included in `data`
    - `map_pips.sh` Script to process all experimental design runs for a single benchmark
    - `process_pips_[A|R].sh` SLURM scripts to process benchmark results for adaptive/randomdesign experiments. Launched by `map_pips.sh` 
    - `process_one.sh` SLURM script to process outputs for a single run
- `scripts/` Benchmarking and post-processing helper scripts
    - `make_sample.R` R script to generate random benchmark datasets. Outputs included in data.Use isolate sheets instead.
    - `make_init_sets.R` R script to generate initial seed batches. Called by `make_*.sh`.
    - `run_[adaptive|random].sh` Replicate runner scripts for adaptive and random designexperiments with varying scale
    - `run_adaptive_old.sh` Replicate runner script for adaptive design experiments with fixedscale
    - `run_rep_[adaptive|random][_old].sh` SLURM launch scripts for replicate runners
    - `run_steps[_adaptive_old].sh` Benchmark launch scripts for dispatching multiple replicaterunners
- `pyseer_comp/` R script and launch scripts to run pyseer on full datasets and to process results
- `out/` Output directory for benchmarks, includes initial seed data.
- `annotate/` Annotation utilities and results processing scripts 
    - `annotate_all.R` R script to annotate all unitigs. Called by `annotate_all_patts.sh`
    - `process_pips_2.R` R script to process results & compute PIPs for a single run. Called byvarious results processing shell scripts
- `plotting/` Shared plotting scripts and `*.gff` reference annotation files
- `analysis_benchmarks/` R scripts to analyze benchmark outputs
    - `scripts/` Scripts to generate manuscript figures and tables
    - `tables/` Output tables
    - `summary_data/` Numeric values for power plots
    - `plots/` Output plots 
    - `mapping/` Annotated GWAS benchmark outputs. Available in zenodo data repository.
- `analysis_lab/` Plotting and analysis scripts for lab application
    - `designs/` Experimental design outputs
    - `mapping/` Annotated GWAS lab outputs. Available in zenodo data repository.
    - `plots/` Plotting output
    - `scripts/` Plotting scripts
    - `run_scripts/` Slurm submission scripts
    

