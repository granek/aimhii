There are several different ways to install AIM-HII:


1.  Docker: recommended for a standard desktop
2.  2: recommended for a bioinformatics server
3.  3: recommended for pros

The bioinformatics software that AIM-HII depends on is described in the 4 section.

# Docker Installation

This is recommended for a standard desktop (Mac or Windows), which is unlikely to have required bioinformatic software installed.  

## Install Docker

1.  See installation instructions for [Mac](https://docs.docker.com/installation/mac/), [Windows](https://docs.docker.com/installation/windows/), [etc](https://docs.docker.com/installation/).
2.  `docker pull granek/aimhii` to download aimhii image
    
    **Note:** If running the following commands on a machine other than Mac or Windows, you may need to use the command `sudo docker` instead of just `docker`.

-   To start Docker daemon from the shell on a Mac

    1.  `boot2docker start` start docker daemon
    2.  ~eval "$(boot2docker shellinit)"~

## Run analysis from manuscript

Both of these commands will output the results to the current directory.  The "subset" results are named `subset_SRR1964709_results.csv`, the "full" results are named `SRR1964709_results.csv`.   An error such as "Are you trying to connect to a TLS-enabled daemon without TLS?", may indicate that docker has not been started.
-   `docker run -v $PWD:/root/aimhii/results -t granek/aimhii make run_subset` to do a test analysis on a subset of the data
-   `docker run -v $PWD:/root/aimhii/results -t granek/aimhii make run_aimhii` to perform the full analysis from the manuscript.

## Analyze your own data

Running `aimhii` from a Docker container requires one step in addition to what you would normally do if it was installed directly on your computer: you have to tell docker where to find the input files, and where to put the results.  If all of the input files are in the current directory, something like the following command will work, saving the results to `results.csv` in the current directory.

    docker run -v $PWD:/mydir \
    -t granek/aimhii aimhii \
    /mydir/genome.fna \
    /mydir/insert.fasta \
    /mydir/adapter.fasta \
    /mydir/R1.fastq.gz \
    /mydir/R2.fastq.gz \
    --outfile /mydir/results.csv

A similar command will work if the input files are in the current directory, or directories within the current directory.  For example if the FASTQ files are in a subdirectory call `fastq_dir`

    docker run -v $PWD:/mydir \
    -t granek/aimhii aimhii \
    /mydir/genome.fna \
    /mydir/insert.fasta \
    /mydir/adapter.fasta \
    /mydir/fastq_dir/R1.fastq.gz \
    /mydir/fastq_dir/R2.fastq.gz \
    --outfile /mydir/results.csv

## Advanced: Access shell in AIMHII docker container

`docker run -i -t granek/aimhii /bin/bash`

# Pip Installation

This is recommended for a server with common bioinformatics software already installed.
This requires that the 4 are already installed.

## Setup a Python Virtual Environment (strongly recommended)

This assumes you have Python [pip](https://pypi.python.org/pypi/pip) and [virtualenv](https://pypi.python.org/pypi/virtualenv) installed.  If [virtualenv](https://pypi.python.org/pypi/virtualenv) is not installed try the command `pip install virtualenv` to install it.  If [pip](https://pypi.python.org/pypi/pip) is not installed, you will need to [install](https://pip.pypa.io/en/stable/installing.html) it first.
1.  `virtualenv aimhii_venv` to set up a virtual environment (see [virtualenv guide](http://docs.python-guide.org/en/latest/dev/virtualenvs/) for details).  **Note:** If the default version of python on your system is older than 2.7, you might need to specify the path to python2.7, for example: `virtualenv -p /usr/local/bin/python2.7 aimhii_venv`
2.  `source aimhii_venv/bin/activate` to enter the virtual environment.  Use the command `deactivate` to leave virtual environment.

## Use pip to install AIMHII

1.  `pip install numpy` (HTSeq needs numpy installed beforehand)
2.  `pip install aimhii`

## To reproduce the analysis from the manuscript (Optional)

1.  `git clone https://granek@bitbucket.org/granek/aimhii.git`
2.  `cd aimhii`
3.  `make run_subset` to do a test analysis on a subset of the data.
4.  `make run_aimhii` to perform the full analysis from the manuscript.

# Clone Repository

This is only recommended for "pros": you are familiar with git and want to explore the source code.  This method requires that the 4 are already installed.

## Download source code

1.  `git clone https://granek@bitbucket.org/granek/aimhii.git`

## To reproduce the analysis from the manuscript (Optional)

1.  `cd aimhii`
2.  `make run_subset` to do a test analysis on a subset of the data.
3.  `make run_aimhii` to perform the full analysis from the manuscript.

# Software Dependencies

The versions given below have been tested with this software, and are known to work, but earlier versions may work perfectly well:

-   [ea-utils](https://code.google.com/p/ea-utils/) (fastq-mcf and fastq-join)
-   [BWA](http://bio-bwa.sourceforge.net/) (version 0.7.5a-r405)
-   [samtools](http://samtools.sourceforge.net/)
-   [python2.7](https://www.python.org/downloads/release/python-279/)
-   [pip](https://pip.pypa.io/en/latest/installing.html) (already included in python version 2.7.9 and higher)

## Only required to replicate analysis from manuscript

-   [gnumake](http://www.gnu.org/software/make/)
-   [git](http://git-scm.com/downloads)
-   [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/#SRA_download.how_do_i_download_and_insta) (version 2.3.5)
-   [curl](http://curl.haxx.se/)

## Python libraries

All of these packages will be installed by `pip aimhii`.  They can be installed separately by running `pip install -r requirements.txt`, using the requirements.txt file included in the repository.

-   biopython
-   HTSeq
-   matplotlib
-   numpy
-   pysam

# Test

[test linking](src/master/docker_install.md)
