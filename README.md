# EXFOR chart

This is an interactive application to inspect the EXFOR database and compare to the predictions of various models. This tool was created as part of the NEUCOS project supported by the European Research Council (ERC - grant No. 646623).

### Installation

1. Download or clone this repository

2. Install the requirements

        python -m pip install pip --upgrade
        python -m pip install -r < requirements.txt

3. Unpack the data file that you received from me to data. The data directory should look like a random collection of files.

4. Run the interactive matplotlib plot:

        python exfor_chart.py -i

### Requirements

- a few GB of hard drive space for `x4i3` package and the data files.
- Python 3 or 2

### Documentation

None. The scientific background behind the choice of models is hidden [in this paper](https://www.nature.com/articles/s41598-017-05120-7).

### Contributions

Any type of contribution or interest is welcome.

### Authors

Anatoli Fedynitch

### License

Code released under [the BSD 3-clause license (see LICENSE)](LICENSE).