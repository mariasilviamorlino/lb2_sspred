Content of the folder
=====================

scripts/
--------

Visualization
* dataviz.py --> functions for plotting, based on matplotlib/seaborn *(yet to write)*

Counters
* taxoparse.py --> it is apart from the others just because I needed specific functions to parse files with taxonomy info
* scop.py --> *DEPENDS ON statstools.py (it imports the compo() function)*
* statstools.py --> a bunch of functions based on dictionary counters. most likely to be imported elsewhere because (at least in my imagination) the functions are as versatile as possible

Parsers
* testsets.py --> to treat the blind test set
* dsspparse.py --> to extract secondary structure string from a dssp output

GOR
* class-gor.py --> implementation of GOR models as a class of their own
* gor-train.py --> train a GOR model from an user-defined training set and store all the parameters (the GOR model) into an output file
* gor-predict.py --> reads GOR parameters from an input model file and predicts the secondary structure for all the sequence in a user-defined testing set, storing the results(predictions) to an output file


### TODO

Organization will be split in packages/modules vs python files actually doing the shit.
The organization will be the following:
* ml
    + class-gor.py
    + dataset_preprocess.py
    + gor-train.py (can be deleted)
    + performance.py
    + sov.py
* stats
    + taxoparse.py
    + scop.py
    + statstools.py
    + testsets.py
+ dataviz.py
+ dsspparse.py

Makefile
--------

### TODO

* Specific makefile for the statistics of the dataset (maybe both blind and training?)
* 
