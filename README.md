This is a set of scripts that I wrote in order to determine some more or less optimal parameters for running KLIP on data from CHARIS. 

The four python scripts are the following:
  - parameterprobing.py: This is the script to run KLIP on a dataset and get contrast and planet detection data.
  - parameter_test_infrastructure.py: This is the majority of the code and the foundation for actualparameterprobing.py.
  - plotting.py: This code is designed to be able to provide a few different types of summary plots which can indicate which parameters are better than others.
  - test_parameter_test_infrastructure.py: This is a testing suite designed to be used with pytest in order to make sure that the majority of the code is bug-free.
  - compilingcontrast.py: a script which can be used for (roughly) finding best contrast curve out of a set of contrast curves
