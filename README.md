# About PETAL (ParallEl paThways AnaLyzer)

PETAL software is written in the Python 3 programming language. It contains a set of tools for pathway analysis and discovery of novel therapeutic targets. The approach allows you to scan and perform a in-depth search of the biological pathway to analyze less recurrent pathways, detect nodes that are far from the initial target nodes and showing the pathway of origin from which it was taken the gene.

![welcome-page](https://github.com/Pex2892/PETAL/blob/master/gui.png)


## Installation

### Download the PETAL repository:

```bash
git clone https://github.com/Pex2892/PETAL.git
```

### Dependencies

PETAL depends on the following libraries:

*   pandas
*   joblib
*   requests
*   beautifulsoup4

A typical user can install the libraries using the following command:

``` bash
pip3 install -r requirements.txt
```

# Testing the Installation

You can test that you have correctly installed the PETAL 
by running the following command:

```bash
python3 main.py
```
After completing an analysis, to view the final output just open the index.html file located in the "demo" folder.


# How to resolve the errors

### bz2 error 
When trying to execute the python script, if the following error occurs:

```bash
from _bz2 import BZ2Compressor, BZ2Decompressor ImportError: No module named '_bz2'
```

You need to install libbz2 and .so files, so that python will be compiled with bz2 support.
```bash
sudo apt install libbz2-dev  # on ubuntu/debian or
sudo yum install bzip2-devel # on fedora

cp /usr/lib/python3.7/lib-dynload/_bz2.cpython-37m-x86_64-linux-gnu.so  /usr/local/lib/python3.7/
```

### lzma error 
When trying to execute the python script, if the following error occurs:

```bash
UserWarning: Could not import the lzma module. Your installed Python is incomplete. 
Attempting to use lzma compression will result in a RuntimeError.
```

You need to install liblzma and .so files, so that python will be compiled with lzma support.
```bash
sudo apt install liblzma-dev  # on ubuntu/debian or
sudo yum install -y xz-devel  # on fedora

cp cp /usr/lib/python3.7/lib-dynload/_lzma.cpython-37m-x86_64-linux-gnu.so /usr/local/lib/python3.7/
```


