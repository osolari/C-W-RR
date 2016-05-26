#! /usr/bin/python

import numpy as np

class lineParser(object):
    def __init__(self, line):
        self.line = line
        columns = line[:-1].split("\t")
        self.ID = "P".join(columns[:2])
        self.info = list()
        #self.gtBinary = self.ID
        self.gq = "" #list()
        self.dp = "" #list()
        self.gtScore = self.ID
        self.maxGQ = 0
        for k in columns[9:]:
            if "1/1:" in k:
                self.info.append(k)
                tmpInfo = k.split(":")
                #self.gq.append(np.float(tmpInfo[1]))
                #self.dp.append(np.float(tmpInfo[2]))
                
                self.gq = self.gq + tmpInfo[1] + "\n"
                self.dp = self.dp + tmpInfo[2] + "\n"
                tmpGQ = np.float(tmpInfo[1])
                if tmpGQ >= self.maxGQ:
                    self.maxGQ = tmpGQ
                if tmpGQ >= 90:
                    self.gtScore = self.gtScore + "\t" + tmpInfo[1]
                    #self.gtBinary = self.gtBinary + "\t" + "1"
                else:
                    self.gtScore = self.gtScore + "\t" + "0"
            else:
                self.gtScore = self.gtScore + "\t" + "0"
                #self.gtBinary = self.gtBinary + "\t" + "0"
        self.gtScore = self.gtScore + "\n"
        #self.gtBinary = self.gtBinary + "\n"  


def vcfHeaderCreator(headerLine):
    columns = headerLine.split("\t")
    row = "P".join(columns[:2]) + "\t" + "\t".join(columns[9:])
    return row[1:]

def vcfLineParser(Line):
    if "1/1:" in Line:
        columns = Line.split("\t")
        row = "P".join(columns[:2])
        GTs = [k.split(":")[0] for k in columns[9:]]
        for j in GTs:
            if j == "1/1":
                row = row + "\t" + "1"
            else:
                row = row + "\t" + "0"
        return row + "\n"
    else:
        return 0

def vcfLineParserINFO(Line):
    if "1/1:" in Line:
        columns = Line[:-1].split("\t")
        row = "P".join(columns[:2])
        
        for k in columns[9:]:
            if "1/1:" in k:
                info
                row = row + "\t" + k
            else:
                row = row + "\t" + "0"
        return row + "\n"
    else:
        return 0


vcf = open("../variants_omid/spatial_variants/dp_4_70_percent_spatial_variants_10000_duplicate_ind_removed_top_293_scaffold_freebayes.recode.vcf","rb")
#vcfParsed = open("../variants_omid/spatial_variants/dp_4_70_percent_spatial_variants_10000_duplicate_ind_removed_top_293_scaffold_freebayes.recode.parsedScores.phred90.tdl.csv", "wb")
#vcfGQ = open("../variants_omid/spatial_variants/dp_4_70_percent_spatial_variants.GQ.tdl.phred90.csv", "wb")
#vcfDP = open("../variants_omid/spatial_variants/dp_4_70_percent_spatial_variants.DP.tdl.phred90.csv", "wb")

phredThreshold = 90

for i in vcf.readlines():
    #if "##" not in i:
    if "#CHROM" in i:
            #vcfParsed.write(vcfHeaderCreator(i))
        header = vcfHeaderCreator(i)
        #elif ("1/1:" in i) and ("0/0" in i):
        #    gtParsed = lineParser(i)
        #    if gtParsed.maxGQ >= phredThreshold:
        #        vcfParsed.write(gtParsed.gtScore)
        #        vcfGQ.write(gtParsed.gq)
        #        vcfDP.write(gtParsed.dp)
vcf.close()
#vcfParsed.close()
#vcfGQ.close()
#vcfDP.close()

