import sys, argparse, os, fnmatch

def main():
    parser = argparse.ArgumentParser(description="This is HTStream python wrapper")
    parser.add_argument("-s", "--sampleSheet", help="This is the sample sheet file.", default="sampleSheet.txt")
    
    args = parser.parse_args()
    parseSampleSheet(args.sampleSheet)

def parseSampleSheet(sampleSheet):
    
    with open(sampleSheet, "r") as ss:
        dataDict = setDataDict()
        apps = []
        sampleDict = {}
        for ssData in ss:
            if ssData[0] == "#":
                #metadata
                nothing = 0
            elif ssData[0] == "*":
                tmpData = ssData[1:].split(":")
                if  tmpData[0] in dataDict:
                    dataDict[tmpData[0]] = os.path.expanduser(tmpData[1].strip())
                elif "app" in tmpData[0]:
                    apps.append(tmpData[1].strip())    
                else:
                    print "ERROR"
                    print ssData
                    print "Not valid * parameter in sample sheet"
                    print "Must be '*appPath', '*rawReads', '*preprocessed', or '*app##' (without quotes and # being literal digits)"
                    exit(1)
            else:
                tmpData = ssData.split("\t")
                if len(tmpData) == 2:
                    sampleDict[tmpData[0]] = tmpData[1].strip()
                else:
                    print "ERROR"
                    print ssData
                    print "Not valid tab delimited samples"
                    print "Please, insure that everything that doesn't have a *, or # in front are tab delimited samples"
                    exit(1)
                #SampleSetup
        findFullAppPath(dataDict["appPath"], apps)
        fastqFiles = createSamplePaths(dataDict["rawReads"], dataDict["preprocessed"], sampleDict)
        setupCommands(apps, fastqFiles)

def setupCommands(apps, fastqFiles):
    for sample in fastqFiles:
        sampleCMD = apps[0]
        #PE1 and PE2
        if sample[0]:
            sampleCMD += " -1 " + ",".join(sample[0]) + " -2 " + ",".join(sample[0]) + " "
        #SE
        if sample[2]:
            sampleCMD += " -U " + ",".join(sample[0]) + " "

        if not sample[0] and not sample[2]:
            print "ERROR"
            print "No FASTQ files found for " + sample[3]
            exit(1)

        sampleCMD += " -O "
        #set up apps 1 through before the last one
        for a in apps[1:-1]:
            sampleCMD += " | " + a + " -O -S "
        
        sampleCMD += " | " + a + " -p " + os.path.join(sample[3], "cleaned_") + " -S "
        print sampleCMD

def getFastqFiles(dirPath, fastqFiles):
    #fastqFiles = [[], [], []]
    for root, dnames, fnames in os.walk(dirPath):
        for fname in fnmatch.filter(fnames, "*.fastq*"):
            if "_R1" in fname:
                tryR2 = fname.replace("_R1", "_R2").strip()
                if os.path.exists(os.path.join(root, tryR2)):
                    fastqFiles[0].append(os.path.join(root, fname))
                    fastqFiles[1].append(os.path.join(root, tryR2))
                else:
                    fastqFiles[2].append(os.path.join(root, fname))
            elif "_R2" in fname:
                nothing = 0
            else:
                fastqFiles[2].append(os.path.join(root, fname))

def createSamplePaths(rawReads, preprocessed, sampleDict):
    #This will have another list appened onto it
    fastqFiles = []
    for k, v in sampleDict.items():
        samples = k.split(',')
        #the list [R1, R2, SE, OUTPUT_PREFIX]
        fastq = [ [], [], [], [] ]
        for sample in samples:
            samplePath = os.path.join(rawReads, sample)
            if os.path.isdir(samplePath):
                getFastqFiles(samplePath, fastq)
            else:
                print "ERROR"
                print samplePath + " : Is not a directory. Make sure you have created your '*rawReads:' directory"

        if not os.path.exists(os.path.join(preprocessed, v)):
            os.makedirs(os.path.join(preprocessed, v))
        fastq[3] = os.path.join(preprocessed, v)
        fastqFiles.append(fastq)
    return fastqFiles
'''
Sets up the applications to be called with their full paths (if need be)
'''
def findFullAppPath(path, apps):
    for root, dnames, fnames in os.walk(path):
        for i in xrange(0, len(apps)):
            a = apps[i].split(' ')
            for fname in fnmatch.filter(fnames, a[0]):
                apps[i] = os.path.join(root, fname)
                if len(a) > 1:
                    apps[i] += " " +  " ".join(a[1:])

def setDataDict():
    dataDict = {}
    dataDict["appPath"] = ""
    dataDict["rawReads"] = ""
    dataDict["preprocessed"] = ""
    return dataDict


if __name__ == "__main__":
    main()
