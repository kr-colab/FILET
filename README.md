This directory contains all of the code necessary to run FILET (Finding Introgressed Loci using Extra-Trees). Briefly, this tool uses an Extra-Trees classifier (Geurts et al. 2005: http://www.montefiore.ulg.ac.be/~ernst/uploads/news/id63/extremely-randomized-trees.pdf) to classify genomic windows from a as one of three different modes of evolution based on population genomic data from a phased two-species sample:
1) Introgression from species 1 to species 2
2) Introgression from species 2 to species 1
3) No introgression

The manuscript describing FILET can be found [here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007341). You can also check out this talk for a brief overview of the method and results from an application to Drosophila simulans and D. sechellia: https://www.youtube.com/watch?v=t6zy8ZsRKt4

In order to run FILET, the user will need several python packages, including scikit-learn (http://scikit-learn.org/stable/), scipy (http://www.scipy.org/), and numpy (http://www.numpy.org/). The simplest way to get all of this up and running is to install Anaconda (https://store.continuum.io/cshop/anaconda/), which contains these and many more python packages and is very easy to install. You will also need GCC (https://gcc.gnu.org/) and GSL (https://www.gnu.org/software/gsl/). The FILET pipeline should then work on any Linux system (and probably OS X too, though this has not yet been tested).

Once you have cloned this repository into the desired directory, cd into this directory and then run:

make all

This will compile the C tools included in this pipeline, and you should then be ready to rock. However, because FILET requires training prior to performing classifications, running the software is a multi-step process. Below, we outline each step, and in examplePipeline.sh we describe the whole process in greater detail, with a complete example run. Here are the steps of the FILET pipeline:

- Simulate training data: Statistics summarizing patterns of variation within each simulated genomic window must be calculated from training data.
- Calculate feature vectors from these training data, and combine into a matrix of labeled training examples.
- Train the FILET classifier.
- Calculate summary statistics from genomic windows we wish to classify.
- Classify genomic data (with known or inferred haplotypic phase).

Not all of these steps are mandatory: if you have generated your own feature vectors from training data, then provided these data are in the proper format (see example in featureVectorsToClassify/) you can proceed to the training step. Similarly, if you have calculated summary statistics from genomic data, then once a classifier has been trained you can skip straight to the data classification step as well (again, provided these statistics are listed in the proper format; see example in dataToClassify/). Note that this pipeline also assumes that you have already simulated training data (with ms-style output), which can be done using msmove (https://github.com/geneva/msmove) or the simulator of your choice. The requirements for simulated data and their file names are described in the comments above Step 4 of examplePipeline.sh.

Note: the full results from our scan in D. simulans and D. sechellia will be added to a directory called simSechResults later today. (I.e. once my dev server comes back up...)
