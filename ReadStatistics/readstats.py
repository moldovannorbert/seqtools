#!/usr/bin/env python3
# Norbert Moldován
import os, re, argparse, pandas as pd, time, statistics, warnings
import seaborn as sns, matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
warnings.simplefilter("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("samPath", help="The path to the directory containing .sam files.")
parser.add_argument("output", help="The path to the outputput folder.")
parser.add_argument("-e","--exclude", dest="exc", type=str, help="Exclude one or more contigs from the analysis. A perfect match of the contig name is required. Contigs must be separated by comas.")
parser.add_argument("-st", "--stat", action="store_true", help="Run only the statistics. For this, specify the path to parsed .txt files in samPath instead of the sam files.")
args = parser.parse_args()

def cigar_interpreter(cigar_string):
    cig_dict = {"M" : 0, "I" : 0, "D" : 0, "N" : 0}
    cigar_list = []
    for number, event, in re.findall("(\d+)([ISMDHN])",cigar_string):
        cigar_list.append((event, number),)
        if event == "M":
            cig_dict["M"] += int(number)
        elif event == "I":
            cig_dict["I"] += int(number)
        elif event == "D":
            cig_dict["D"] += int(number)
        elif event == "N":
            cig_dict["N"] += int(number)
    return cigar_list, cig_dict

def statCalc(outPath, filename):
    df = pd.read_csv("".join((outPath,filename)), header=None, sep="\t")
    df = df.sort_values([3], ascending=False).drop_duplicates(subset=[0]).sort_index()
    try:
        excL = (args.exc).split(",")
        excDf = df[~df[1].isin(excL)]
    except:
        excDf = df
    excDf = excDf.sort_values([3], ascending=False).drop_duplicates(subset=[0]).sort_index()
    
    readCount = len(df[0])
    mappedCount = len(excDf[excDf[3] > 0])
    mappedFrac = mappedCount/readCount
    avgLength = df[2].mean()
    
    fracMismatch = ((excDf[excDf[3] > 0][4]).sum()/(excDf[excDf[3] > 0][3]).sum())
    fracI = ((excDf[excDf[3] > 0][5]).sum()/(excDf[excDf[3] > 0][3]).sum())
    fracD = ((excDf[excDf[3] > 0][6]).sum()/(excDf[excDf[3] > 0][3]).sum())
    fracGC = ((excDf[excDf[3] > 0][7]).sum()/(excDf[excDf[3] > 0][3]).sum())
    try:
        stDevLength = statistics.stdev(df[2])
        stDevMismatch = statistics.stdev((excDf[excDf[3] > 0][4]))
        stDevI = statistics.stdev((excDf[excDf[3] > 0][5]))
        stDevD = statistics.stdev((excDf[excDf[3] > 0][6]))
        stGC = statistics.stdev((excDf[excDf[3] > 0][7]))
    except:
        print("Only one datapoint.")
        stDevLength = 0
        stDevMismatch = 0
        stDevI = 0
        stDevD = 0
        stGC = 0
        
    avgLengthMapped = (excDf[excDf[3] > 0][3]).mean()
    try:
        stDevLengthMapped = statistics.stdev(excDf[excDf[3] > 0][3])
    except:
        print("Only one datapoint.")
        stDevLengthMapped = 0
        
    medLength = df[2].median()
        
    medLengthMapped = (excDf[excDf[3] > 0][3]).median()
    maxLength = df[2].max()
    maxLengthMapped = (excDf[excDf[3] > 0][3]).max()
    minLength = df[2].min()
    minLengthMapped = (excDf[excDf[3] > 0][3]).min()
    
    with open("".join((args.output, "ReadStat.tsv")), 'a') as outputputfile:
        outputputfile.write("\t".join((filename.strip(".sam.txt"), str(readCount), 
                                       str(mappedCount), str(mappedFrac), str(minLength), 
                                       str(minLengthMapped), str(avgLength), str(stDevLength),
                                       str(avgLengthMapped), str(stDevLengthMapped), str(medLength),
                                       str(medLengthMapped), str(maxLength), str(maxLengthMapped), 
                                       str(fracMismatch), str(stDevMismatch), str(fracI), str(stDevI),
                                       str(fracD), str(stDevD), str(fracGC), str(stGC), "\n")))
    outputputfile.close()

def loader(filesInPath, txtPath, filenameInPath):
    print("Progressing to: ",filenameInPath)
    with open("".join((txtPath, filenameInPath, ".txt")), 'a') as outputputfile:    
        with open("".join((filesInPath,filenameInPath))) as sam:
            for line in sam:
                if line and not line.startswith("@"):
                    line = line.split("\t")
                    readID = line[0]
                    contig = line[2]
                    cigar_list, cig_dict = cigar_interpreter(line[5])
                    read_length = len(line[9])
                    if contig != "*":
                        map_length = int(cig_dict["M"])+int(cig_dict["D"])-1
                        NM_tag = line[11]
                        ins_count = int(cig_dict["I"])
                        del_count = int(cig_dict["D"])
                        try:
                            startSoft = int(line[5].split("S")[0])
                        except:
                            startSoft = 0
                        AlignedSeq = line[9][startSoft:startSoft+int(cig_dict["M"])+ins_count]
                        GC_count = int(AlignedSeq.count("G") + AlignedSeq.count("C"))
                        GC_frac = GC_count/(int(cig_dict["M"])+int(cig_dict["I"]))
                        mismatch_count = int(NM_tag.split(":")[2]) - (ins_count + del_count)
                    else:
                        map_length = 0
                        mismatch_count = 0
                        ins_count = 0
                        del_count = 0
                        GC_count = 0
                        GC_frac = 0
                        AlignedSeq = "*"
                    outputputfile.write("\t".join((readID, str(contig), str(read_length), str(map_length), str(mismatch_count), str(ins_count), str(del_count), str(GC_count), str(GC_frac), "\n")))
    outputputfile.close()
    
def Plotter(dataPath, OutPath):
    PlotDf = pd.DataFrame(columns = ["ReadID", "Contig", "Read Length (nt)", "Mapped Length (nt)", "Mismatch count", "Insertion count", "Deletion count", "GC count", "GC frac", "Sample"])
    
    for files in os.walk(dataPath):
        for filename in files[2]:
            if filename.endswith(".txt"):
                HDf = pd.DataFrame()
                HDf = pd.read_csv("".join((dataPath, filename)), header=None, sep="\t")
                HDf = HDf.sort_values([3], ascending=False).drop_duplicates(subset=[0]).sort_index()
                HDf.columns = ["ReadID", "Contig", "Read Length (nt)", "Mapped Length (nt)", "Mismatch count", "Insertion count", "Deletion count", "GC count", "GC frac", "*"]
                HDf.drop(columns = "*", inplace = True)
                HDf["Sample"] = filename.strip(".sam.txt")
            try:
                PlotDf = PlotDf.append(HDf)
            except:
                pass
            
    ReadStDf = pd.read_csv("".join((OutPath, "ReadStat.tsv")), index_col = 0, sep="\t")
    
    MapDf = PlotDf[PlotDf["Mapped Length (nt)"] != 0]
    MapDf["Mapped Length (nt)"] = MapDf["Mapped Length (nt)"].astype('float64')
    MapDf["Read Length (nt)"] = MapDf["Read Length (nt)"].astype('float64')
    MapDf["Mismatch frac"] = (MapDf["Mismatch count"]/MapDf["Mapped Length (nt)"]).astype('float64')
    MapDf["Insertion frac"] = (MapDf["Insertion count"]/MapDf["Mapped Length (nt)"]).astype('float64')
    MapDf["Deletion frac"] = (MapDf["Deletion count"]/MapDf["Mapped Length (nt)"]).astype('float64')
    
    LenDist = sns.kdeplot(PlotDf["Read Length (nt)"], shade = True)
    LenDist = sns.kdeplot(MapDf["Mapped Length (nt)"], shade = True)
    sns.despine()
    LenDist.set(title = "Mapped/Unmapped Read Length Distribution", xlabel='Read length', ylabel='Distribution',
                xscale = "log")
    LenDist.figure.savefig("".join((OutPath, "/figs/MappedUnmappedLength.svg")))
    plt.clf()
    
    sns.set(style="white", palette="muted", color_codes=True)
    
### Aligned read length plotted on one plot
    for s in set(MapDf.Sample):
        MapAlignedDist = sns.distplot(MapDf[MapDf.Sample == s]["Mapped Length (nt)"], label = s, hist = False)
    sns.despine()
    plt.legend()
    MapAlignedDist.set(title = "Aligned Read Length Distribution", xlabel='Read length', ylabel='Distribution',
                       xscale = "log")
    MapAlignedDist.figure.savefig("".join((OutPath, "/figs/AlignedReadLength_dist.svg")))
    plt.clf()
### Total read length plotted on one plot
    for s in set(MapDf.Sample):
        MapTotalDist = sns.distplot(MapDf[MapDf.Sample == s]["Read Length (nt)"], label = s, hist = False)
    sns.despine()
    plt.legend()
    MapTotalDist.set(title = "Raw Read Length Distribution", xlabel='Read length', ylabel='Distribution',
                     xscale = "log")
    MapTotalDist.figure.savefig("".join((OutPath, "/figs/RawReadLength_dist.svg")))
    plt.clf()
### GC% distribution - Short
    for s in set(MapDf.Sample):
        GC_plot = sns.distplot(MapDf[MapDf.Sample == s]["GC frac"], label = s, hist = False)
    sns.despine()
    plt.legend()
    GC_plot.set(title = "GC distribution", xlabel='GC fraction', ylabel='Distribution')
    GC_plot.figure.savefig("".join((OutPath, "/figs/GC_dist.svg")))
    plt.clf()
    
### Aligned read length violinplot
    MapAlignedDist=sns.violinplot(x = MapDf["Sample"], y = MapDf["Mapped Length (nt)"], dataset = MapDf)
    sns.set(style="white", font_scale=0.8)
    sns.despine()
    MapAlignedDist.set_xticklabels(MapAlignedDist.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    MapAlignedDist.figure.savefig("".join((OutPath, "/figs/AlignedReadLength_violin.svg")))
    plt.clf()
### Total read length violinplot      
    MapTotalDist=sns.violinplot(x = MapDf["Sample"], y = MapDf["Read Length (nt)"], dataset = MapDf)
    sns.set(style="white", font_scale=0.8)
    sns.despine()
    MapTotalDist.set_xticklabels(MapTotalDist.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    MapTotalDist.figure.savefig("".join((OutPath, "/figs/RawReadLength_violin.svg")))
    plt.clf()
### GC% distribution - Long
    GC_plot=sns.violinplot(x = MapDf["Sample"], y = MapDf["GC frac"], dataset = MapDf)
    sns.set(style="white", font_scale=0.8)
    sns.despine()
    GC_plot.set_xticklabels(GC_plot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    GC_plot.figure.savefig("".join((OutPath, "/figs/GC_dist_violin.svg")))
    plt.clf()
### Mismatch distribution - violin
    MmPlot=sns.violinplot(x = MapDf["Sample"], y = MapDf["Mismatch frac"], dataset = MapDf, linewidth=0.3)
    sns.set(style="white", font_scale=0.8)
    sns.despine()
    MmPlot.set_xticklabels(MmPlot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    MmPlot.figure.savefig("".join((OutPath, "/figs/Mismatch_violin.svg")))
    plt.clf()
### Insertion distribution - violin
    InsPlot=sns.violinplot(x = MapDf["Sample"], y = MapDf["Insertion frac"], dataset = MapDf, linewidth=0.3)
    sns.set(style="white", font_scale=0.8)
    sns.despine()
    InsPlot.set_xticklabels(InsPlot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    InsPlot.figure.savefig("".join((OutPath, "/figs/Insertion_violin.svg")))
    plt.clf()
### Deletion distribution - violin
    DelPlot=sns.violinplot(x = MapDf["Sample"], y = MapDf["Deletion frac"], dataset = MapDf, linewidth=0.3)
    sns.set(style="white", font_scale=0.8)
    sns.despine()
    DelPlot.set_xticklabels(DelPlot.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.tight_layout()
    DelPlot.figure.savefig("".join((OutPath, "/figs/Deletion_violin.svg")))
    plt.clf()

### Pie chart of total read counts  
    SeqPie = ReadStDf.plot.pie(y= "Read count", autopct='%1.1f%%', pctdistance=0.8, figsize=(7, 5), labels=None)
    SeqPie.legend(labels=ReadStDf.index, loc="upper right", bbox_to_anchor=(1.3, 1))
    SeqPie.set(title = "Total Read Count", ylabel="")
    SeqPie.figure.savefig("".join((args.output, "/figs/TotalReadCount_pie.svg")))
    plt.clf()
### Pie charts of mapped read counts
    SeqPie = ReadStDf.plot.pie(y= "Mapped read count", autopct='%1.1f%%', pctdistance=0.8, figsize=(7, 5), labels=None)
    SeqPie.legend(labels=ReadStDf.index, loc="upper right", bbox_to_anchor=(1.3, 1))
    SeqPie.set(title = "Mapped Read Count", ylabel="")
    SeqPie.figure.savefig("".join((OutPath, "/figs/MappedReadCount_pie.svg")))
    plt.clf()

## Mapped read length and mismatch count per read scatterplot - short
    # Scatter = sns.scatterplot(data = MapDf, x="Mapped Length (nt)", y="Mismatch count", hue="Sample", linewidth=0, marker=".", alpha=0.3, fit_reg=True)
    # sns.set(style="white", font_scale=0.5)
    # sns.despine()
    # Scatter.figure.savefig("".join((args.output, "/figs/ReadLength_mmatch_short.svg")))
    # plt.clf()

#### Mapped read length and mismatch count per read scatterplot - long
#    Scatter = sns.FacetGrid(MapDf, row = "Sample", height=1.2, aspect=2)
#    Scatter.map(sns.scatterplot,  "Mapped Length (nt)", "Mismatch count", linewidth=0)
#    sns.set(style="white", font_scale=0.8)
#    sns.despine()
#    plt.subplots_adjust(top=0.9)
#    Scatter.fig.suptitle("Mapped read length and mismatch count per read")
#    Scatter.savefig("".join((OutPath, "/figs/RLength_mmatch_long.svg")))
#    plt.clf()



def Main():
    Tt = time.time()
    processes = []
    samFileNames = []
    txtFileNames = []
    pool = Pool(processes = os.cpu_count()-1)
            
    if not args.stat:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
            os.makedirs("".join((args.output,"/figs")))
            os.makedirs("".join((args.output,"/data")))
        print("\n\nCollecting variables from SAM files.\n")
        
        for files in os.walk(args.samPath):
            for filename in files[2]:
                samFound = False
                if filename.endswith(".sam"):
                    samFound = True
                    samFileNames.append(filename)
    
        if samFound == False:
            print("\n\n No SAM files found in ", args.samPath,"\n")
            pass
            
        t=time.time()
        pool.map(partial(loader, args.samPath, "".join((args.output,"/data/"))), samFileNames)
        pool.close()
        pool.join()
        print("Collecting variables finished in", time.time()-t,"s")
    
        with open("".join((args.output, "ReadStat.tsv")), 'w') as outputputfile:
            outputputfile.write("\t".join(("Samples", "Read count", "Mapped read count", 
                                           "Mapped fraction", "Min read length", "Min mapped read length", 
                                           "Mean read length", "Stdev length", "Mean mapped read length", 
                                           "Stdev mapped length", "Median read length", "Median mapped read length", 
                                           "Max read length", "Max mapped read length", "Mismatch fraction", "Stdev mismatch",
                                           "Insertion fraction", "Stdev insertion", "Deletion fraction", "Stdev deletion",
                                           "GC fraction of aligned", "Stdev GC content", "\n")))
        outputputfile.close()
        
        pool = Pool(processes = os.cpu_count()-1)
        
        print("\n Calculating statistics.")
        
        t=time.time()
        
        for files in os.walk("".join((args.output,"/data/"))):
            for filename in files[2]:
                fileFound = 0
                if filename.endswith(".txt"):
                    fileFound += 1
                    txtFileNames.append(filename)
       
        pool.map(partial(statCalc, "".join((args.output,"/data/"))), txtFileNames)
        pool.close()
        pool.join()
        
        print("\n Calculating statistics finished in", time.time()-t,"s")
        
        print("\n Preparing plots.")
        t=time.time()
        Plotter("".join((args.output,"data/")), args.output)
        print("\n Plotting finished in", time.time()-t,"s")
        
    else:
        if not os.path.exists(args.output):
            os.makedirs(args.output)
            os.makedirs("".join((args.output,"/figs")))

        with open("".join((args.output, "ReadStat.tsv")), 'w') as outputputfile:
            outputputfile.write("\t".join(("Samples", "Read count", "Mapped read count", 
                                           "Mapped fraction", "Min read length", "Min mapped read length", 
                                           "Mean read length", "Stdev length", "Mean mapped read length", 
                                           "Stdev mapped length", "Median read length", "Median mapped read length", 
                                           "Max read length", "Max mapped read length", "Mismatch fraction", "Stdev mismatch",
                                           "Insertion fraction", "Stdev insertion", "Deletion fraction", "Stdev deletion",
                                           "GC fraction of aligned", "Stdev GC content", "\n")))
        outputputfile.close()
        
        pool = Pool(processes = os.cpu_count()-1)
        
        print("\n Calculating statistics.")
        
        t=time.time()
        
        for files in os.walk(args.samPath):
            for filename in files[2]:
                fileFound = 0
                if filename.endswith(".txt"):
                    fileFound += 1
                    txtFileNames.append(filename)
       
        pool.map(partial(statCalc, args.samPath), txtFileNames)
        pool.close()
        pool.join()
        
        print("\n Calculating statistics finished in", time.time()-t,"s")
        
        print("\n Preparing plots.")
        t=time.time()
        Plotter(args.samPath, args.output)
        print("\n Plotting finished in", time.time()-t,"s")
        print("\n Total run took:", time.time()-Tt,"s") 
    
    
if __name__ == "__main__":
	Main()	