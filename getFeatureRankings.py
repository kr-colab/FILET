import sys, os
import numpy as np
import random
from sklearn.externals import joblib
from sklearn.ensemble import ExtraTreesClassifier

classifierPickleFileName = sys.argv[1]

statsToUse, header, grid_search = joblib.load(classifierPickleFileName)

importances = grid_search.best_estimator_.feature_importances_

# Print the feature ranking
print("Feature ranking from an ExtraTreesClassifier:")

outLines = []
for f in range(len(importances)):
    outLines.append((importances[f], "feature %s; value: %g" % (header[f+1], importances[f])))
outLines.sort()
outLines.reverse()
featureRank = 1
for importance, outLine in outLines:
    print("%d. " %(featureRank) + outLine)
    featureRank += 1
print("done\n")
