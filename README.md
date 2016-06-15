<center>
# CoLoRMap
</center>

## Installation
In order to install CoLoRMap, you should first fetch the source code from CoLoRMap git repository.

```bash
git clone --recursive https://github.com/sfu-compbio/colormap.git
```

After obtaining the code, you need to install the dependencies. CoLoRMap uses *BWA*, *SAMtools*, and *Minia*. In order to build these dependencies, change to the source directory `colormap` and use `make deps` command.

```bash
cd colormap
make deps
```

At last, you can compile CoLoRMap binaries simply by running `make` command.

```bash
make
```

## Correcting long reads
To correct long reads, you can use **runCorr.sh** script:

```bash
./runCorr.sh <pacbio.fasta> <illumina.fastq> <outPrefix> <threads>
```
After finishing this, the corrected long reads are stored in `<outPrefix>_corr.fasta` file in `<outPrefix>` directory.

## Improving the correction using One-End Anchors (OEAs) ##
The script **runOEA.sh** can be used to further improve the quality of corrected long reads by using One-End Anchors (OEAs) to extend the borders of the corrected regions.
```bash
./runOEA.sh <pacbio_corr.fasta> <illumina.fastq> <outPrefix> <threads>
```
When this is done, the corrected long reads are stored in `<outPrefix>_oea.fasta` file in `<outPrefix>` directory.

## Pulication
CoLoRMap: Correcting Long Reads by Mapping short reads<br/>
*Haghshenas E, Hach F, Sahinalp SC, Chauve C*<br/>
Accepted in ECCB 2016, to appear in Bioinformatics

## Contact
Please report problems and bugs on [issues page](https://github.com/sfu-compbio/colormap/issues). Otherwise, contact `ehaghshe[at]sfu[dot]ca`
