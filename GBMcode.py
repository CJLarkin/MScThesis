#main data analysis

#\t for each cell
#\n for new line

import csv
from itertools import combinations
import numpy as np
import scipy
from scipy.stats.stats import pearsonr #returns Pearson's correlation coefficient, 2-tailed p-value
from scipy.stats.stats import spearmanr #returns Spearman's correlation coefficient, 2-tailed p-value
import seaborn as sns
import operator
import matplotlib.pyplot as plt
import matplotlib
import math

#returns mutation, gene, study, sample
def read(file, empty):
    with open(file, newline = '') as GBM_data:
        lineID = 0
        headings = []
        mutations = []
        lines = csv.reader(GBM_data, delimiter='\n')
        for line in lines:
            cell = csv.reader(line, delimiter='\t')
            for info in cell:
                foo = info[0]
                bar = info[1]
                columnID = 0
                if lineID == 0:
                    headings.append(info)
                    titles = headings[0]
                    title1 = titles[0]
                    title2 = titles[1]
                for data in info:
                    output = []
                    if (data != empty) & (lineID != 0) & (columnID != 0) & (columnID != 1):
                        #print ("Mutation: ", data,"Gene: ", titles[columnID], title1 ,": ", foo, title2,": ", bar)
                        output.append(data)
                        output.append(titles[columnID])
                        output.append(foo)
                        output.append(bar)
                        #print (output)
                        mutations.append(output)
                    columnID += 1
            lineID +=1
    return mutations

#outputs similar genes
def getGenes(data, file):
    genes = []
    for i in data:
        gene = i[1]
        genes.append(gene)
    genes = list(dict.fromkeys(genes))
    
    #print (genes)
    with open(file, newline = '') as RNA_data:
        geneID = []
        lines = csv.reader(RNA_data, delimiter='\n')
        for line in lines:
            cell = csv.reader(line, delimiter='\t')
            for info in cell:
                columnID = 0
                for data in info: #going through cell by cell
                    if columnID == 0:
                        for name in genes:
                            if name in data:
                                geneID.append(data)
                columnID += 1
    return geneID

#returns [[gene, expression level 1, 2, 3,...], [gene2, exp1, 2, 3,...],...]
def RNAexp(genes, file):
    output = []
    with open(file, newline = '') as RNA_data:
        lines = csv.reader(RNA_data, delimiter='\n')
        for line in lines:
            cell = list(csv.reader(line, delimiter='\t'))
            for info in cell:
                for i in genes:
                    if info[0] == i:
                        things = []
                        for data in info:
                            things.append(data)
                        output.append(things)
    return output

def writeFile(file, data):
    f = open(file,'w')
    f.write("Gene" + "\t" + "Log2 Expression Level")
    f.write("\n")
    for i in data:
        count = 0
        for n in i:
            if count != 0:
                try:
                    f.write("%s\t" % float(n))
                except ValueError:
                    f.write("%s\t" % n)
            else:
                f.write("%s\t" % n)
            count +=1
        f.write("\n")
    f.close()
    
def toFloat(data): #takes in singular set of [gene, exp1, exp2,...], returns ['gene', [exp1,exp2,...]] with exp as float
    output = []
    gene = data[0]
    for i in data[1:]:
        try:
            output.append(float(i))
        except ValueError:
            output.append(0) #this might screw with the results, let's see...
            gene = "null"
    return gene, output

#takes in output of RNAexp, ie [[gene, expression level 1, 2, 3,...], [gene2, exp1, 2, 3,...],...]
def correlation(data):
    comb = list(combinations(np.arange(0,len(data),1),2))
    relationship = {}
    for n in comb:
        gene1, exp1 = toFloat(data[n[0]])
        gene2, exp2 = toFloat(data[n[1]])
        corr, p = scipy.stats.pearsonr(exp1, exp2)
        relationship[gene1 + ", " + gene2] = corr
    return relationship

def corrRanked(data):
    comb = list(combinations(np.arange(0,len(data),1),2))
    relationship = {}
    for n in comb:
        gene1, exp1 = toFloat(data[n[0]])
        gene2, exp2 = toFloat(data[n[1]])
        corr, p = scipy.stats.spearmanr(exp1, exp2)
        relationship[gene1 + ", " + gene2] = corr
    return relationship

#writes GINML file with specified nodes and edges
def WriteGINML(corr_data, outfile):
        PathwayID = "MAPK_corr" #just name of system
        ProteinList = nodes(corr_data)
        AllEdges = edges(corr_data)

        f = open('{0}.ginml'.format(outfile), 'w')
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE gxl SYSTEM "http://ginsim.org/GINML_2_2.dtd">\n')
        f.write('<gxl xmlns:xlink="http://www.w3.org/1999/xlink">\n')
        lineadd = '  <graph class="regulatory" id="' + PathwayID + '" nodeorder="' + ' '.join(ProteinList) + '">\n'
        f.write(lineadd)
        f.write('\n')
        f.write('    <nodestyle background="#ffffff" foreground="#000000" text="#000000" shape="RECTANGLE" width="45" height="25"/>\n')
        f.write('    <nodestyle background="#ffffff" foreground="#000000" text="#000000" shape="RECTANGLE" width="45" height="25"/>\n')
        f.write('    <edgestyle color="#000000" pattern="SIMPLE" line_width="1" properties="positive:#00c800 negative:#c80000 dual:#0000c8"/>\n')
        f.write('    <edgestyle color="#000000" pattern="SIMPLE" line_width="1" properties="positive:#00c800 negative:#c80000 dual:#0000c8"/>\n')
        f.write('\n')

        for i in range(len(ProteinList)):
            lineadd2 = '    <node id="' + ProteinList[i] + '" maxvalue="1">\n'
            f.write(lineadd2)
            f.write('      <value val="1">\n' + '      </value>\n')

            r = 300                                     #300 is an arbitrary value for now, will make it a function of node number
            x = r + (-1)*r*math.cos((i/len(ProteinList))*2*math.pi)   #This defines position of the x coordinate
            y = r + (-1)*r*math.sin((i/len(ProteinList))*2*math.pi)   #This defines position of the y coordinate
            f.write('      <nodevisualsetting x="' + str(round(x)) + '"' + ' y="' + str(round(y)) + '" style=""/>\n') #Make a circular spread using the sin(360/len(ProteinList), make x y variable
            f.write('    </node>\n')                                           #Turn x and y into variables

        f.write('\n')

        for i in range(len(AllEdges)):
            lineadd3 = '    <edge id="' + AllEdges[i][0] + ':' + AllEdges[i][1] + '" from="' + AllEdges[i][0] + '" to="' + AllEdges[i][1] + '" minvalue="1" sign="' + AllEdges[i][2] + '">\n'
            f.write(lineadd3)
            f.write('      <edgevisualsetting style=""/>\n')
            f.write('    </edge>\n')

        f.write('\n')
        f.write('  </graph>\n')
        f.write('</gxl>')
        f.close()

        
#takes in result of corr functions and returns list of genes invovled in relationships with r >= 0.6
def nodes(data):
    nodes = []
    for i in range(len(data)):
        #print (abs(data[i][1]))
        if 1 > abs(data[i][1]) >= 0.6:
            split = [x.strip() for x in sorted_corrsR[i][0].split(',')]
            #print (split[0])
            if (split[0] != "null") & (split[1] != "null"):
                nodes.append(split[0])
                nodes.append(split[1])
    nodes = list(dict.fromkeys(nodes))
    return nodes

#takes in result of corr functions and returns list of lists of gene1, gene2, positive/negative
def edges(data):
    edges = []
    for i in range(len(data)):
        edge = []
        if 1 > data[i][1] >= 0.6:
            split = [x.strip() for x in sorted_corrsR[i][0].split(',')]
            if (split[0] != "null") & (split[1] != "null"):
                edge.append(split[0])
                edge.append(split[1])
                edge.append("positive")
                edges.append(edge)
        if -1 < data[i][1] <= -0.6:
            split = [x.strip() for x in sorted_corrsR[i][0].split(',')]
            if (split[0] != "null") & (split[1] != "null"):
                edge.append(split[0])
                edge.append(split[1])
                edge.append("negative")
                edges.append(edge)
    return edges


#calling the functions
data = read("/Users/collinjlarkin/MSc_project/GenomeData/mutations.txt", "NA")
file = "/Users/collinjlarkin/MSc_project/GenomeData/GBMnormalized.txt"
destFile = "/Users/collinjlarkin/MSc_project/GenomeData/test.txt"
commonGenes = getGenes(data, file)
rna = RNAexp(commonGenes, file)
#print (rna[0])
writeFile(destFile,rna)


corrsR = corrRanked(rna)
sorted_corrsR = sorted(corrsR.items(), key=operator.itemgetter(1))

nodes(sorted_corrsR)
edges(sorted_corrsR)

WriteGINML(sorted_corrsR,"MAPKtest")

#frequency plots
count = 0
for i in range(len(rna)):
    gene1, exp1 = toFloat(rna[i])
    if gene1 != "null":
        count += 1
        e1 = sns.kdeplot(exp1, shade="True")
plt.title("Gene Expression")
plt.xlabel("Normalized Gene Expression Level")
plt.ylabel("Frequency")
plt.show()
print (count)

avg = []
for i in range(len(rna)):
    gene1, exp1 = toFloat(rna[i])
    if gene1 != "null":
        avg.append(np.mean(exp1))
e1 = sns.kdeplot(avg, shade="True")
plt.title("Gene Expression")
plt.xlabel("Normalized Gene Expression Level")
plt.ylabel("Frequency")
plt.show()

#low genes
lowGenes = []
for i in range(len(rna)):
    gene, exp = toFloat(rna[i])
    if gene != 'null':
        for n in exp:
            if n < 5.5:
                lowGenes.append(gene)
    lowGenes = list(dict.fromkeys(lowGenes))

print (lowGenes)

#high genes
highGenes = []
for i in range(len(rna)):
    gene, exp = toFloat(rna[i])
    if gene != 'null':
        for n in exp:
            if n > 9.5:
                highGenes.append(gene)
    highGenes = list(dict.fromkeys(highGenes))

print (highGenes)

#bimodal genes
def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3

intersection(lowGenes,highGenes)



