# module load Python/2.7.15-foss-2018b
import sys
import os
import gzip
import numpy as np

# this script filter SComatic output based on the control somatic mutations from WGS data

sId = sys.argv[1]
tId = sys.argv[2]
print "Control ID:",tId

# somatic
print("Control: somatic")
control="/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mutationCalling/hg38BulkRes/snvs_ICGC_%s_somatic_snvs_conf_8_to_10.hg38.bed" % tId

assert os.path.exists(control)
assert "somatic" in control or "germline" in control

# snATAC
inFile = "/b06x-isilon/b06x-m/mbCSF/results/humanTumor/mbSpatial/atacMut/%s/Step4_VariantCalling/%s.calling.step2.tsv.gz" % (sId,sId)

print "Input:",inFile

targCellTypes =  ["C1","C2","C3"]
print("Targ cell types:",targCellTypes)
 

refSet = set()

for line in open(control):
    if line.startswith("#"):
        continue
    items = line.strip().split("\t")
    snvId = "%s:%s" % (items[0].replace("chr",""),items[1])   
    refSet.add(snvId)

print "Selection from control bulk:",len(refSet)
#print refSet


assert ".tsv.gz" in inFile
resFile = open(inFile.replace(".tsv.gz",".somatic.tsv"), "w" )

novel = 0
found = 0
vafStats = {}
targCellTypeIds = [] # only focus on target cell types
 
for line in gzip.open(inFile):

    items = line.strip().split("\t")
 
    if line.startswith("#CHROM"):
        line = line.replace("\n", "\tmeanVAF\n"  )
        resFile.write(line)
        for i in range(len(items)):
            if items[i] in targCellTypes:
                targCellTypeIds.append(i)
        print "Targ cell type IDs:", targCellTypeIds        
        continue

    if line.startswith("##"):
        continue
    snvId = "%s:%s" % (items[0].replace("chr",""),items[1])   
    novel += 1
    #print snvId   

    #if snvId in refSet: # somatic
    if not (snvId in refSet): # non-germline
        #print "Found:", snvId
        if not items[0].startswith("chr"):
            continue
        
        fInfo = items[5]
        cType = items[6]
        if ("Multi-allelic" in fInfo) or ("|" in items[4]): # several somatic changes
            continue
            
        # skip unclear specificity
        if (cType == ".")  or ("Normal" in cType):
            continue         
        
        targElems = [items[i] for i in targCellTypeIds]
        # compute VAF
        str_array = items[14].split(",")
        str_float = [float(x) for x in str_array]
        vaf = np.mean(  str_float   )
        #print("%s-%s %s %.4f" % (items[0], items[1],cType, vaf)  )
        if cType in vafStats:
            vafStats[cType].append( vaf )
        else:
            vafStats[cType] = [ vaf ]
        found += 1
        line = line.replace("\n", "\t%.4f\n" % vaf )

        resFile.write(line)

print "Found in target total: %d, passed control bulk: %d" % (novel,found) 
print "VAF stats:"
for cT in vafStats:
    vafs = vafStats[cT]    
    print cT, len(vafs),  np.mean(vafs)









