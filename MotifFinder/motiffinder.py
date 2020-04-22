#!/usr/bin/env python3
# Norbert Moldován
# import warnings
# warnings.simplefilter("ignore")
import os.path, pandas as pd, argparse, time
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from functools import partial
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument("--genome_path", "-g", dest="genome_path", type=str, required=True, help="The path to the genome .fasta file.")
parser.add_argument("--annotation_path", "-i", dest="annotation_path", type=str, required=True, help="The path to the transcript annotation directory containing .gff/.gff3 files.")
parser.add_argument("--out_path", "-s", dest="out_path", type=str, required=True, help="The path to the output folder.")
parser.add_argument("--tss_int", "-ts", default=5, type=int, help="The interval from wich the sequence of TSS should be called.")
parser.add_argument("--tes_int", "-te", default=50, type=int, help="The distance from the TES where the GU richness should be analyzed")
parser.add_argument("--min_readcount", "-mrc", default=10, type=int, help="The minimal read-count necesarry to write a feature in the TSS/TES.fa files.")
parser.add_argument("--exclude", "-e", default="", type=str, help="Exclude a contig.")
args = parser.parse_args()

colNames = ["contig",
           "source",
           "feature",
           "start",
           "end",
           "readCount",
           "strand",
           "cds",
           "name",
           "gc_dist",
           "gc_seq",
           "ca_dist",
           "ca_seq",
           "ta_dist",
           "ta_seq",
           "pa_dist",
           "pa_seq",
           "tes_seq",
           "tss_seq",
           "TSS_A_frac",
           "TSS_C_frac",
           "TSS_G_frac",
           "TSS_T_frac",
           "TES_A_frac",
           "TES_C_frac",
           "TES_G_frac",
           "TES_T_frac",
           "samples"]
### Checks if output folder exists. If not, creates it. 
def DirManager():
    if not os.path.exists(args.out_path):
        os.makedirs(args.out_path)

### Loads sequence from fasta file.
def SeqLoader(record, posL, posR, strand):
    sequence = ""
    if strand == "-":
        #print("SeqLoad",posL,posR)
        pos = SeqFeature(FeatureLocation((posL-1), posR), type="gene")
        sequence = str(record[pos.location.start:pos.location.end].reverse_complement().seq)

    elif strand == "+":
        pos = SeqFeature(FeatureLocation(posL, posR), type="gene")
        sequence = str(record[pos.location.start:pos.location.end].seq)
    #print(posL,posR,sequence)
    return sequence    
   
###The Finder function searches for the motifs (motif) in the genomic intervall defined by the position of the TSS or TES (row) and the ideal intervall of the motif (p1, p2),
###and creates a .gff output.
def Finder(AnnotDf, motiffile, cont):
    seq = str()
    annotDf = AnnotDf.loc[AnnotDf.contig == cont].copy()
    annotDf["ca_seq"] = annotDf["ca_seq"].astype(str)
    annotDf["gc_seq"] = annotDf["gc_seq"].astype(str)
    annotDf["ta_seq"] = annotDf["ta_seq"].astype(str)
    annotDf["pa_seq"] = annotDf["pa_seq"].astype(str)
    annotDf["tss_seq"] = annotDf["tss_seq"].astype(str)
    annotDf["tes_seq"] = annotDf["tes_seq"].astype(str)    
    #print(annotDf[["contig","start","end"]].head())
    for record in SeqIO.parse(args.genome_path, "fasta"):
        if record.id == cont:
            for ai, ar in annotDf.iterrows():
                if ar["strand"] == "-":
                    if ar["feature"] in ["mRNA","tss","TSS"]:
                        for mi, mr in motiffile[motiffile.type == "GC_box"].iterrows():
                            gcS = ar["end"] + mr["start"]
                            gcE = ar["end"] + mr["end"]
                            seq = SeqLoader(record, gcS, gcE, "-")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].find(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"gc_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"gc_seq"] = mr["motif"]
                        
                        for mi, mr in motiffile[motiffile.type == "CCAAT_box"].iterrows():
                            caS = ar["end"] + mr["start"]
                            caE = ar["end"] + mr["end"]
                            seq = SeqLoader(record, caS, caE, "-")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].find(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"ca_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"ca_seq"] = mr["motif"]
                        
                        for mi, mr in motiffile[motiffile.type == "TATA_box"].iterrows():
                            taS = ar["end"] + mr["start"]
                            taE = ar["end"] + mr["end"]
                            seq = SeqLoader(record, taS, taE, "-")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].find(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"ta_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"ta_seq"] = mr["motif"]
                    
                        tsS = ar["end"] - args.tss_int
                        tsE = ar["end"] + args.tss_int
                        annotDf.at[ai,"tss_seq"] = tssSeq = SeqLoader(record, tsS, tsE, "-")
                        
                        cA = (tssSeq.upper()).count("A")
                        cG = (tssSeq.upper()).count("G")
                        cC = (tssSeq.upper()).count("C")
                        cT = (tssSeq.upper()).count("T")
                        
                        if len(tssSeq) > 0:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = cA/len(tssSeq)
                            annotDf.at[ai,"".join(("TES_G_frac"))] = cG/len(tssSeq)
                            annotDf.at[ai,"".join(("TES_C_frac"))] = cC/len(tssSeq)
                            annotDf.at[ai,"".join(("TES_T_frac"))] = cT/len(tssSeq)
                        else:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_G_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_C_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_T_frac"))] = "nan"
                    
                    if ar["feature"] in ["mRNA","tes","TES"]:
                        for mi, mr in motiffile[motiffile.type == "PAS"].iterrows():
                            paS = ar["start"] + mr["start"]
                            paE = ar["start"] + mr["end"]
                            seq = SeqLoader(record, paS, paE, "-")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].find(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"pa_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"pa_seq"] = mr["motif"]
                    
                        teS = ar["start"] - args.tes_int
                        teE = ar["start"] + args.tes_int
                        annotDf.at[ai,"tes_seq"] = tesSeq = SeqLoader(record, teS, teE, "-")
                        
                        cA = (tesSeq.upper()).count("A")
                        cG = (tesSeq.upper()).count("G")
                        cC = (tesSeq.upper()).count("C")
                        cT = (tesSeq.upper()).count("T")
                        
                        if len(tesSeq) > 0:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = cA/len(tesSeq)
                            annotDf.at[ai,"".join(("TES_G_frac"))] = cG/len(tesSeq)
                            annotDf.at[ai,"".join(("TES_C_frac"))] = cC/len(tesSeq)
                            annotDf.at[ai,"".join(("TES_T_frac"))] = cT/len(tesSeq)
                        else:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_G_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_C_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_T_frac"))] = "nan"
                else:
                    if ar["feature"] in ["mRNA","tss","TSS"]:
                        for mi, mr in motiffile[motiffile.type == "GC_box"].iterrows():
                            gcE = ar["start"] - mr["start"]
                            gcS = ar["start"] - mr["end"]
                            seq = SeqLoader(record, gcS, gcE, "+")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].rfind(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"gc_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"gc_seq"] = mr["motif"]
                        
                        for mi, mr in motiffile[motiffile.type == "CCAAT_box"].iterrows():
                            caE = ar["start"] - mr["start"]
                            caS = ar["start"] - mr["end"]
                            seq = SeqLoader(record, caS, caE, "+")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].rfind(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"ca_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"ca_seq"] = mr["motif"]
                        
                        for mi, mr in motiffile[motiffile.type == "TATA_box"].iterrows():
                            taE = ar["start"] - mr["start"]
                            taS = ar["start"] - mr["end"]
                            seq = SeqLoader(record, taS, taE, "+")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].rfind(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"ta_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"ta_seq"] = mr["motif"]
                                
                        tsS = ar["start"] - (args.tss_int+1)
                        tsE = ar["start"] + args.tss_int
                        annotDf.at[ai,"tss_seq"] = tssSeq = SeqLoader(record, tsS, tsE, "+")
                        
                        cA = (tssSeq.upper()).count("A")
                        cG = (tssSeq.upper()).count("G")
                        cC = (tssSeq.upper()).count("C")
                        cT = (tssSeq.upper()).count("T")
                        
                        if len(tssSeq) > 0:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = cA/len(tssSeq)
                            annotDf.at[ai,"".join(("TES_G_frac"))] = cG/len(tssSeq)
                            annotDf.at[ai,"".join(("TES_C_frac"))] = cC/len(tssSeq)
                            annotDf.at[ai,"".join(("TES_T_frac"))] = cT/len(tssSeq)
                        else:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_G_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_C_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_T_frac"))] = "nan"
                        
                    if ar["feature"] in ["mRNA","tes","TES"]:
                        for mi, mr in motiffile[motiffile.type == "PAS"].iterrows():
                            paE = ar["end"] - mr["start"]
                            paS = ar["end"] - mr["end"]
                            seq = SeqLoader(record, paS, paE, "+")
                            motifStart = seq[0:mr["end"] - len(mr["motif"])].rfind(mr["motif"])
                            if motifStart > -1:
                                annotDf.at[ai,"pa_dist"] = mr["end"]-motifStart
                                annotDf.at[ai,"pa_seq"] = mr["motif"]
                    
                        teS = ar["end"] - (args.tes_int+1)
                        teE = ar["end"] + args.tes_int
                        annotDf.at[ai,"tes_seq"] = tesSeq = SeqLoader(record, teS, teE, "+")

                        cA = (tesSeq.upper()).count("A")
                        cG = (tesSeq.upper()).count("G")
                        cC = (tesSeq.upper()).count("C")
                        cT = (tesSeq.upper()).count("T")
                        
                        if len(tesSeq) > 0:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = cA/len(tesSeq)
                            annotDf.at[ai,"".join(("TES_G_frac"))] = cG/len(tesSeq)
                            annotDf.at[ai,"".join(("TES_C_frac"))] = cC/len(tesSeq)
                            annotDf.at[ai,"".join(("TES_T_frac"))] = cT/len(tesSeq)
                        else:
                            annotDf.at[ai,"".join(("TES_A_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_G_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_C_frac"))] = "nan"
                            annotDf.at[ai,"".join(("TES_T_frac"))] = "nan"
    
    with open("".join((args.out_path,"Data.tsv")), 'a') as outputfile:
        annotDf.to_csv(outputfile, index=False, header=colNames, sep="\t")
    outputfile.close()

def FileLoader():
    print("Loading annotation files... \n")
    nTrDf = pd.DataFrame(columns=colNames)
    helperDf = pd.DataFrame(columns=colNames)
    
    for files in os.walk(args.annotation_path):  
        for filename in files[2]:
            if filename.endswith(".gff"):
                clean_filename=filename.split(".gff")
                helperDf = pd.read_csv("".join((files[0], clean_filename[0],".gff")), names=colNames, sep = "\t", comment="#") #Append each file after another
                helperDf["samples"] = helperDf["samples"].astype("str")
                for hI, hR in helperDf.iterrows():
                    helperDf.at[hI,"samples"] = clean_filename[0]
                nTrDf = nTrDf.append(helperDf, ignore_index = True, sort=False)
            elif filename.endswith(".gff3"):
                clean_filename=filename.split(".gff3")
                helperDf = pd.read_csv("".join((files[0], clean_filename[0],".gff3")), names=colNames, sep = "\t", comment="#") #Append each file after another
                helperDf["samples"] = helperDf["samples"].astype("str")
                for hI, hR in helperDf.iterrows():
                    helperDf.at[hI,"samples"] = clean_filename[0]
                nTrDf = nTrDf.append(helperDf, ignore_index = True, sort=False)
    
    nTrDf = nTrDf.sort_values(by=["contig","start","strand"])
    nTrDf = nTrDf.reset_index(drop = True)
    #print(nTrDf)
    return nTrDf

def Main():
    tS = time.time()
    print("\n")
    print("=====Welcome to motif finder=====")
    motiffile = pd.read_csv("motif.txt", names=["motif","type","start","end"] ,header=None, comment='#', sep="\t")

    annotDf = FileLoader()
    annotDf = annotDf.dropna(how = "all").sort_values(by=["contig","start","end","strand"]).reset_index(drop = True)
    
    DirManager()
    
    contigSet = set(annotDf["contig"])
    
    if args.exclude in contigSet:
        contigSet.remove(args.exclude)  
        
    open("".join((args.out_path,"Data.tsv")), 'w').close()
    
    processes = []
    pool = Pool(processes = os.cpu_count()-1)
    func = partial(Finder, annotDf, motiffile)
    pool.map(func, contigSet)
    pool.close()
    pool.join()
    tE = time.time()
    
    annotDf = pd.read_csv("".join((args.out_path,"Data.tsv")), names=colNames, sep = "\t", comment="#")
    annotDf.drop_duplicates(colNames,keep = "first",inplace = True)
    annotDf = annotDf.sort_values(by=["contig","feature","start","strand"])
    print("The run took",tE-tS)

if __name__ == "__main__":
	Main()	