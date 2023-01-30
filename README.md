# UbiqPred
Our project goal was to develop a model that finds Ubiquitin binding domain (UBDs).
We refined a an end-to-end, interpretable geometric deep learning model called ScanNet https://github.com/jertubiana/ScanNet
In this repository there are a couple of modules:
dataCreation.py-preprocessing of the data set.
cath.py - Dividing the preprocessed data(chains) to 5 homolous sublists.
histogram.py - creating an histogram for the UBDs probabillities and non UBDs probabillities
transfer_learning_train.py - a module from scanNet with a couple of changes.
UBDs - directory that contains the .cif files of our data set.
datasets/BCE-our 5 sublists for the cross validation.
