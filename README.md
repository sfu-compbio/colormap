# CoLoRMap: Correcting Long noisy Reads by Mapping short reads

## Installation
In order to install CoLoRMap, you should first fetch the source code from CoLoRMap git repository.

```bash
git clone --recursive https://github.com/sfu-compbio/colormap.git
```

Please note that command argument `--recursive` is necessary for downloading submodules automatically. After obtaining the code, you need to install the dependencies. CoLoRMap uses *BWA*, *SAMtools*, and *Minia*. In order to build these dependencies, change to the source directory `colormap` and use `make deps` command.

```bash
cd colormap
make deps
```

At last, you can compile CoLoRMap binaries simply by running `make` command.

```bash
make
```

CoLoRMap corrects long reads in two different steps: (i) using a shortest path (SP) algorithm. (ii) using an One-End Anchor (OEA) based algorithm. 

## Preparing the short read data
SP algorithm does not need paired-end information, but OEA algorithm actually uses paired-end information. In both cases, the program expects to be fed with a single short read file. In case of OEA algorithm, program expects paired-end short reads in interleaved/interlaced format. 

Usually, paired-end short reads are stored in two different files. A single interleaved/interlaced read file can be obtained using `fastUtils` program which can be found in `bin` directory after building the program:

```bash
cd testData
../bin/fastUtils shuffle -1 ill_1.fastq -2 ill_2.fastq -o ill.fastq
```

## Correcting long reads
To correct long reads, you can use **runCorr.sh** script:

```bash
../runCorr.sh pac.fasta ill.fastq testCorr pre 4
```
This runs shortest path correction algorithm for long reads stored in `pac.fasta` by short reads stored in `ill.fastq` using 4 threads. When this is done, the corrected long reads are stored in `testCorr/pre_sp.fasta` file.

## Improving the correction using One-End Anchors (OEAs) ##
The script **runOEA.sh** can be used to further improve the quality of corrected long reads by using One-End Anchors (OEAs) to extend the borders of the corrected regions.
```bash
../runOEA.sh testCorr/pre_sp.fasta ill.fastq testOEA pre 4
```
This runs OEA algorithm for pre-corrected long reads stored in `testCorr/pre_sp.fasta` by paired-end short reads stored in interleaved/interlaced format in `ill.fastq` using 4 threads. When this is done, the corrected long reads are stored in `testOEA/pre_oea.fasta` file.

## Publication
*Haghshenas E., Hach F., Sahinalp S.C. and Chauve C.*, "CoLoRMap: Correcting Long Reads by Mapping short reads" [Bioinformatics](http://bioinformatics.oxfordjournals.org/content/32/17/i545.short) (2016) 32 (17): i545-i551
DOI: [10.1093/bioinformatics/btw463](http://dx.doi.org/10.1093/bioinformatics/btw463)

## Contact
Please report problems and bugs on [issues page](https://github.com/sfu-compbio/colormap/issues). Otherwise, contact `ehaghshe[at]sfu[dot]ca`
