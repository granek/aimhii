1.  Download this repository 
    -   Using git: `git clone https://bitbucket.org/granek/aimhii.git`
    -   Alternatively, download a zip file: <https://bitbucket.org/granek/aimhii/downloads>
2.  Install the prerequisites detailed below.
3.  Move into the repository directory: `cd aimhii`
4.  Run the analysis: `make`

# Prerequisites

## Software

The versions given below have been tested with this repository, and are known to work, but earlier versions may work perfectly well
-   [BWA](http://bio-bwa.sourceforge.net/) (version 2.2.3)
-   [NCBI SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) (version 2.3.5)

-   Likely to be installed on most Linux servers

    -   [Python 2.7](https://www.python.org/download/releases/2.7.8/)
    -   [curl](http://curl.haxx.se/download.html)

-   Python libraries

    All of these packages can be installed using [pip](https://pypi.python.org/pypi/pip) by running `pip install -r requirements.txt`, using the requirements.txt file included in the repository.
    -   biopython
    -   HTSeq
    -   matplotlib
    -   numpy
    -   pysam
    
    You might consider installing and running in a [Virtual Environment](https://pypi.python.org/pypi/virtualenv) (see [virtualenv guide](http://docs.python-guide.org/en/latest/dev/virtualenvs/)), in which case the following commands will get you running :
    1.  `git clone https://bitbucket.org/granek/aimhii.git`
    2.  `virtualenv venv`
        -   **Note:** If the default version of python on your system is older than 2.7, you might need to specify the path to python2.7, for example: `virtualenv -p /usr/local/bin/python2.7 venv`
    3.  `source venv/bin/activate`
    4.  `cd aimhii`
    5.  `pip install -r info/requirements.txt`
    6.  `make`

## Raw Sequence data

Running `make` will automatically download the [[<http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56353][raw> data] from this work is available from [GEO](http://www.ncbi.nlm.nih.gov/geo/).
