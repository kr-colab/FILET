import os, random, sys

trainingOutDir = sys.argv[1]
os.system("mkdir -p %s" %(trainingOutDir))

"""
AIC: 10303.1237186
with u = 3.500000e-09
Nref: 487835.088398
nu1_0 : 117605.344694
nu2_0 : 4878347.23386
nu1 : 9279970.1758
nu2 : 26691.779717
T : 86404.6219829
2Nref_m12 : 0.0128753943002
2Nref_m21 : 0.0861669095413
"""

sampleSize1 = int(sys.argv[3])
sampleSize2 = int(sys.argv[4])
numSites=10000
thetaMean, thetaOverRhoMean, nu1Mean, nu2Mean, m12Times2Mean, m21Times2Mean, TMean = 68.29691232, 0.2, 19.022761, 0.054715, 0.025751, 0.172334, 0.664194


#msmove 34 1000 -t 68.29691232 -r 341.4845616 10000 -I 2 20 14 -n 1 19.022761 -n 2 0.054715 -eg 0 1 6.576808 -eg 0 2 -7.841388 -ma x 0.025751 0.172334 x -ej 0.664194 2 1 -en 0.664194 1 1
trainingSampleNumber = int(sys.argv[2])

def drawUnif(m, fold=0.5):
    x = m*fold
    return random.uniform(m-x, m+x)

def drawParams(thetaMean, thetaOverRhoMean, nu1Mean, nu2Mean, m12Times2Mean, m21Times2Mean, TMean, mig=False):
    theta = drawUnif(thetaMean)
    thetaOverRho = drawUnif(thetaOverRhoMean)
    rho = theta/thetaOverRho
    nu1 = drawUnif(nu1Mean)
    nu2 = drawUnif(nu2Mean)
    T = drawUnif(TMean)
    m12 = drawUnif(m12Times2Mean)
    m21 = drawUnif(m21Times2Mean)
    if mig:
        migTime = random.uniform(0, T/4)
        migProb = 1-random.random()
        return theta, rho, nu1, nu2, T, m12, m21, migTime, migProb
    else:
        return theta, rho, nu1, nu2, T, m12, m21

def writeTbsFile(params, outFileName):
    with open(outFileName, "w") as outFile:
        for paramVec in params:
            outFile.write(" ".join([str(x) for x in paramVec]) + "\n")
            
if __name__ == '__main__':
    slurm_cmd = 'sbatch -t 08:00:00 --output=/dev/null --error=/dev/null --mem=8G --wrap "{0} | tee {1} && gzip {1}"'

    for simType, sampleNumber, outDir in [("train", trainingSampleNumber, trainingOutDir)]:
        noMigParams, mig12Params, mig21Params = [], [], []
        for i in range(sampleNumber):
            theta, rho, nu1, nu2, splitTime, m12, m21 = drawParams(thetaMean, thetaOverRhoMean, nu1Mean, nu2Mean, m12Times2Mean, m21Times2Mean, TMean)
            #paramVec = [theta, rho, nu1, nu2, m12, m21, splitTime, splitTime]
            paramVec = [theta, rho, nu1, nu2, 0, 0, splitTime, splitTime]
            noMigParams.append(paramVec)
    
            theta, rho, nu1, nu2, splitTime, m12, m21, migTime, migProb = drawParams(thetaMean, thetaOverRhoMean, nu1Mean, nu2Mean, m12Times2Mean, m21Times2Mean, TMean, mig=True)
            paramVec = [theta, rho, nu1, nu2, 0, 0, splitTime, splitTime, migTime, migProb]
            mig12Params.append(paramVec)
    
            theta, rho, nu1, nu2, splitTime, m12, m21, migTime, migProb = drawParams(thetaMean, thetaOverRhoMean, nu1Mean, nu2Mean, m12Times2Mean, m21Times2Mean, TMean, mig=True)
            paramVec = [theta, rho, nu1, nu2, 0, 0, splitTime, splitTime, migTime, migProb]
            mig21Params.append(paramVec)
            
        noMigTbsFileName = "%s/noMig.tbs" %(outDir)
        noMigSimCmd = "msmove/gccRelease/msmove %d %d -t tbs -r tbs %d -I 2 %d %d -n 1 tbs -n 2 tbs -eg 0 1 6.576808 -eg 0 2 -7.841388 -ma x tbs tbs x -ej tbs 2 1 -en tbs 1 1 < %s" %(sampleSize1+sampleSize2, sampleNumber, numSites, sampleSize1, sampleSize2, noMigTbsFileName)
        writeTbsFile(noMigParams, noMigTbsFileName)
        os.system(slurm_cmd.format(noMigSimCmd, os.path.join(outDir, 'noMig.msOut')))
        #os.system("echo \"%s > %s/noMig.msOut\" | qsub -q HP8g7,HPblg7,galaxy -w e -N msmove -o /dev/null -e /dev/null" %(noMigSimCmd, outDir))
    
        mig12TbsFileName = "%s/mig12.tbs" %(outDir)
        mig12SimCmd = "msmove/gccRelease/msmove %d %d -t tbs -r tbs %d -I 2 %d %d -n 1 tbs -n 2 tbs -eg 0 1 6.576808 -eg 0 2 -7.841388 -ma x tbs tbs x -ej tbs 2 1 -en tbs 1 1 -ev tbs 1 2 tbs < %s" %(sampleSize1+sampleSize2, sampleNumber, numSites, sampleSize1, sampleSize2, mig12TbsFileName)
        writeTbsFile(mig12Params, mig12TbsFileName)
        os.system(slurm_cmd.format(mig12SimCmd, os.path.join(outDir, 'mig12.msOut')))
    
        mig21TbsFileName = "%s/mig21.tbs" %(outDir)
        mig21SimCmd = "msmove/gccRelease/msmove %d %d -t tbs -r tbs %d -I 2 %d %d -n 1 tbs -n 2 tbs -eg 0 1 6.576808 -eg 0 2 -7.841388 -ma x tbs tbs x -ej tbs 2 1 -en tbs 1 1 -ev tbs 2 1 tbs < %s" %(sampleSize1+sampleSize2, sampleNumber, numSites, sampleSize1, sampleSize2, mig21TbsFileName)
        writeTbsFile(mig21Params, mig21TbsFileName)
        os.system(slurm_cmd.format(mig21SimCmd, os.path.join(outDir, 'mig21.msOut')))
