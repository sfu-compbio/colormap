# CoLoRMap

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
path-to-CoLoRMap/bin/fastUtils shuffle -1 read_1.fastq -2 read_2.fastq > read_paired.fastq
```

## Correcting long reads
To correct long reads, you can use **runCorr.sh** script:

```bash
./runCorr.sh <pacbio.fasta> <illumina.fastq> <outDirectory> <outPrefix> <threads>
```
After finishing this, the corrected long reads are stored in `<outPrefix>_sp.fasta` file in `<outDirectory>` directory.

## Improving the correction using One-End Anchors (OEAs) ##
The script **runOEA.sh** can be used to further improve the quality of corrected long reads by using One-End Anchors (OEAs) to extend the borders of the corrected regions.
```bash
./runOEA.sh <pacbio_corr.fasta> <illumina.fastq> <outDirectory> <outPrefix> <threads>
```
When this is done, the corrected long reads are stored in `<outPrefix>_oea.fasta` file in `<outDirectory>` directory.

## Publication
CoLoRMap: Correcting Long Reads by Mapping short reads<br/>
*Haghshenas E, Hach F, Sahinalp SC, Chauve C*<br/>
Bioinformatics 32.17 (2016): i545-i551

## Contact
Please report problems and bugs on [issues page](https://github.com/sfu-compbio/colormap/issues). Otherwise, contact `ehaghshe[at]sfu[dot]ca`
