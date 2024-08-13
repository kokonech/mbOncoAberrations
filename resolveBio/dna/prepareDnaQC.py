import os.path
import sys
import fnmatch
import glob

def find_fq_files(directory_path, target_extension):
    found_files = []

    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith(target_extension):
                found_files.append(os.path.join(root, file))

    return found_files

# this script generates code for submitting jobs per cell to server

sId = sys.argv[1] # sID "MB248"
inPath = sys.argv[2] # path to fastq files

# depends on OS
resDir = "/omics/odcf/analysis/OE0290_projects/bioskryp/ResolveOme/Medulloblastoma/dna/Result/BJ-DNA-QC/%s" % sId

scriptName ="submit_DNA_QC_%s.sh" % sId
scriptFile = open(scriptName, "w")

scriptFile.write("SID=%s\n\n" % sId)


directories = [d for d in os.listdir(inPath) if os.path.isdir(os.path.join(inPath, d))]

sInfo = {} 

for cId in directories:

    if cId == "0_all":
        continue

    print(cId)

  
    cPath = find_fq_files("%s/%s" % (inPath,cId) , "_R1.fastq.gz")
    assert(len(cPath) == 1)
    
    f1Path = cPath[0]
    if  not os.path.exists(f1Path):
        print("Error! File not found:", f1Path)
        sys.exit(-1)

    # original DNA QC
    scriptFile.write("CID=%s\n" % (cId)  )
    scriptFile.write("FILEPATH=%s\n" % (f1Path.replace("_R1.fastq.gz", "") ) )
    scriptFile.write("qsub -N BJ_DNA_QC_${SID}_${CID} -v TUMOR=\"$SID\",INPUT=\"$CID\",READS=\"$FILEPATH\"  run_BJ_DNA_QC.sh\n\n")
 
