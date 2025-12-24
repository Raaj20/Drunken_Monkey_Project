# Drunken_Monkey_Project

# Install Homebrew 
% git clone https://github.com/Homebrew/brew homebrew
% eval "$(homebrew/bin/brew shellenv)"
brew update --force --quiet
chmod -R go-w "$(brew --prefix)/share/zsh"

# Create an environment in Conda
% brew install miniforge
% conda update -n base -c conda-forge conda

% conda --version
conda 25.11.1

conda create -n mirna_pipeline python=3.10

% echo $SHELL
/bin/zsh

conda init zsh

# Restart Terminal

% conda activate mirna_pipeline

# Install required packages
conda install -c bioconda fastqc
conda install -c bioconda bbmap
conda install -c bioconda bowtie
conda install -c bioconda samtools

 % fastqc --version                
FastQC v0.12.1

BBDuk --version
BBTools version 39.52

% bowtie --version
/Users/viraa/homebrew/Caskroom/miniforge/base/envs/mirna_pipeline/bin/bowtie-align-s version 1.3.1
64-bit
Built on Mac-1733801069421.local
2024-12-10T04:50:10

% samtools --version
samtools 1.23
Using htslib 1.23

# Indexing and setting up Directories
head mature.fa 
grep '^>' mature.fa | head -n 20

## Check if Macaca mullata is abbreviated as mml
grep '^>' mature.fa | sed 's/>//' | cut -d'-' -f1 | sort -u

### Double-checking if there are any Rhesus miRNAs to align
grep -i 'Macaca mulatta' mature.fa | head
grep '^>' maturemml.fa | head -n 50

awk '/^>/{p = ($0 ~ /^>mml/)} p{print}' mature.fa > mature_mml.fa

grep -c '^>' mature.fa 
48885

grep -c '^>' mature_mml.fa
912

grep '^>' mature_mml.fa | head -n 50
>mml-miR-200c-5p MIMAT0026560 Macaca mulatta miR-200c-5p
>mml-miR-200c-3p MIMAT0002195 Macaca mulatta miR-200c-3p

### Convert Us to Ts and doublecheck number of miRNAs before and after conversion
awk '/^>/{print; next} { gsub(/[Uu]/, "T"); print }' mature_mml.fa > mature_mmlTs.fa

grep -c '^>' mature_mmlTs.fa 
912

head mature_mmlTs.fa                                                                
>mml-miR-200c-5p MIMAT0026560 Macaca mulatta miR-200c-5p
CGTCTTACCCAGCAGTGTTTGG
>mml-miR-200c-3p MIMAT0002195 Macaca mulatta miR-200c-3p
AATACTGCCGGGTAATGATGGA


### Indexing the Rhesus macaque + U-->T converted .fa files
bowtie-build mature_mmlTs.fa mature_mmlTs
Settings:
  Output files: "mature_mmlTs.*.ebwt"


