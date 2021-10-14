#!/usr/bin/python
#=================================================================================================================================
# Name                  : fastQTLresult.addinfo.GTEx_TSS_pLI.v1.py
# Created On            : 2021-2-2
# Author                : lili
# Comment               :
# Version               : 0.1
# Last Modified By      : lili
# Last Modified On      : 2021.2.2 add the information of GTEx,TSS,pLI for eQTL
# Modified              :
# Description           :
#=================================================================================================================================
def safe_open(file,mode):
    file = os.path.abspath(file)
    if not os.path.exists(file):
        exit('%s is not exists.' %file)
    elif file.endswith('.gz'):
        import gzip
        return gzip.open(file,mode)
    else:
        return open(file,mode)

def getpos(name,title,name2=''):
    ntitle = [i.lower() for i in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    elif name2 != '' and name2.lower() in ntitle:
        pos = ntitle.index(name2.lower())
    else:
        if name2 != '':
            exit('Warning: %s and %s not in title.' %(name,name2))
        else:
            exit('Warning: %s not in title.' %name)
    return pos

#get information from GTEx for eSNP_eGene
#chr_pos_ref_alt_GeneID Pvalue_effect
def getInfo(infile):
    print("dicInfo start: ")
    inf = infile.split(",")[0]
    tag = infile.split(",")[1]
    dicInfo = {}
    for ii in safe_open(inf,'r'):
        if ii.startswith("variant_id"):
            iiline = ii.strip().split("\t")
            pid = getpos('variant_id',iiline)
            pgene = getpos('gene_id',iiline)
            ppvalue = getpos('pval_nominal',iiline)
            peffect = getpos('slope',iiline)
            dicInfo['title'] = '\t'.join([tag+"_pvalue",tag+"_effect"])
        else:
            iitmp = ii.strip().split('\t')
            sid = iitmp[pid].split("_b37")[0]
            gene = iitmp[pgene].split(".")[0]
            marker = '_'.join([sid,gene])
            pvalue = iitmp[ppvalue]
            effect = iitmp[peffect]
            if not dicInfo.get(marker,False):
                dicInfo[marker] = "\t".join([pvalue,effect])
            else:
	            print('This marker has more than one pvalue! '+marker)
    print("dicInfo done!")
    return dicInfo

#get info from TSS (transcrip start pos from refGene) for eGene_TSS?
#genename	chromsome	TranscripStart
#WASH7P	1	14361
#-->{genename:chr_TSS}
def getTSS(infile):
    print("get the TSS info of the genename from refGene")
    dicInfo = {}
    with safe_open(infile,'r') as inf:
        for i in inf:
            if i.startswith("genename"):
                dicInfo["title"] = "_".join(["chr","TSS"])
            else:
                ii = i.strip().split("\t")
                gene = ii[0]
                chrom = ii[1]
                tss = ii[2]
                if not dicInfo.get(gene,False):
                    dicInfo[gene] = "_".join([chrom,tss])
                else:
                    if "_" in chrom:continue
                    elif chrom == dicInfo[gene].split("\t")[0] and float(tss)<float(dicInfo[gene].split("_")[1]):
                        dicInfo[gene] = "_".join([chrom,tss])
                        print("The gene has more than one tss, please check! The old tss is "+dicInfo[gene]+"; the new tss is "+"_".join([chrom,tss]))
                    else:
                        continue
        print("Done! TSS info!")
        return dicInfo

#get pLI from gnomAD v2
#gene	pLI	gene_id	chromosome	start_position	end_position
#MED13	1.0000e+00	ENSG00000108510	17	60019966	60142643
#-->{genename:pLI}
def getpLI(infile):
    print("get the pLI info from gnomAD V2")
    dicInfo = {}
    with safe_open(infile,'r') as inf:
        for i in inf:
            if i.startswith("gene"):
                dicInfo["title"] = "pLI"
            else:
                ii = i.strip().split("\t")
                gene = ii[0]
                pLI = ii[1]
                if not dicInfo.get(gene,False):
                    dicInfo[gene] = pLI
                else:
                    print('This gene has more than one pLI! '+ gene)
        print("Done! pLI info!")
        return dicInfo

#get the eGene genename from the FPKM data geneinfo
#geneID	Gene name	Interpro Description	chromosome	start	end	length	gene_biotype
#ENSG00000252498	--	--	6	76353800	76353893	94	snRNA
#ENSG00000175826	CTDNEP1	Dullard phosphatase domain, eukaryotic||HAD-like domain||-||FCP1 homology domain	1771469107155810	998	protein_coding
##-->{geneID:genename}
def geteGeneName(infile):
    print("get the eGene's genename from FPKM data geneinfo")
    dicInfo = {}
    with safe_open(infile,'r') as inf:
        for i in inf:
            if i.startswith("geneID") or i.startswith("gene_id"):
                dicInfo["title"] = "eGeneName"
            else:
                ii = i.strip().split("\t")
                gene = ii[0]
                genename = ii[1]
                if not dicInfo.get(gene,False):
                    dicInfo[gene] = genename
                else:
                    print('This geneID has more than one genename! '+ gene)
        print("Done! eGene genename info!")
        return dicInfo

#get genename and TSS for eGene
def geteGeneName_TSS(geneinfofile,dictss):
    print("get the eGene genename from FPKM data geneinfo and the TSS of eGENE")
    with safe_open(geneinfofile,'r') as geneinf1:
        for i in geneinf1:
            if i.startswith("geneID") or i.startswith("gene_id"):
                dicInfo["title"] = "eGeneName\tchr_tss"
            else:
                ii = i.strip().split("\t")
                gene = ii[0]
                genename = ii[1]
                #for tss info
                info = []
                info.append(genename)
                if genename == "--":
                    info.append("._.")
                elif dictss.get(genename,False):
                    info.append(dictss[genename])
                else:
                    print("This genename has no TSS info: "+genename)
                    info.append("._.")

                if not dicInfo.get(gene,False):
                    dicInfo[gene] = "\t".join(info)
                else:
                    print('This geneID has more than one genename! '+ gene)
        print("Done! genename and TSS info of eGene!")
        return dicInfo


#add info
def addInfo(infile,dicInfo1,dicInfo2,dicgenename,dicpLI,dicTSS,outfile):
#add the Info to annovar result
#fastQTL result header:
# gene_id	gene_name	gene_chr	gene_start	gene_end	strand	num_var	beta_shape1	beta_shape2	true_df	pval_true_df	variant_id	tss_distancechr	variant_pos	ref	alt	num_alt_per_site	rs_id_dbSNP	ma_samples	ma_count	maf	ref_factor	pval_nominal	slope	slope_se	pval_perm	pval_beta	qval	pval_nominal_threshold
    with safe_open(infile,'r') as inf, open(outfile,'w') as out:
        print("Add the information of Info start!")
        for line in inf:
            if line.startswith("SNP"):
      	        title = line.strip().split("\t")
                psnp = getpos('SNP',title)#eSNP
   	            pgene = getpos('gene',title)#eGene
                pFDR = getpos('FDR',title)
                pgenename = getpos('GeneName', title)#eSNP_genename

                #eSNP,eGene,eGeneName,Beta,pvalue,FDR,GTExinfo,eGenepLI,eGeneTSS,CHR,POS,REF,ALT,ID,QUAL,FILTER ,eSNP_GeneName,eSNP_genepLI,eSNP_geneTSS,Description,Func,...
	        out.write("\t".join(["eSNP","eGene","eGene_GeneName","Beta","pvalue","FDR"] + dicInfo1["title"].split("\t")+ dicInfo2["title"].split("\t")+ ["eGene_pLI","eGene_TSS","length_POS_eGeneTSS","eSNP_CHR","eSNP_POS","eSNP_REF","eSNP_ALT","ID","QUAL","FILTER","eSNP_GeneName","eSNP_genepLI","eSNP_geneTSS","length_POS_eSNPGeneTSS"]+ title[(pgenename+1)::])+"\n")
            else:
                snp_line = line.strip().split("\t")
                esnp = snp_line[psnp]
                egene = snp_line[pgene]
                esnp_GeneName = snp_line[pgenename]
                outsnpline = []

                #eSNP-eGene info
                egene_GeneName = dicgenename[egene]
                outsnpline.extend([esnp,egene,egene_GeneName]+snp_line[pgene+1:pFDR+1])

                ##add GTEx information -Sigmoid
 	        marker = "_".join([esnp,egene])
	        GTEx_InfoList = ["."]*(len(dicInfo1["title"].split("\t")))
                if dicInfo1.get(marker,False):
                    GTEx_InfoList[0:len(GTEx_InfoList)] = dicInfo1[marker].split("\t")
                else:
                    GTEx_InfoList[0:len(GTEx_InfoList)] = ["."]*2
                    print("The variant-gene has no Info: "+marker)
                ##add GTEx information2 -Transverse:
                GTEx_InfoList2 = ["."]*(len(dicInfo2["title"].split("\t")))
                if dicInfo2.get(marker,False):
                    GTEx_InfoList2[0:len(GTEx_InfoList2)] = dicInfo2[marker].split("\t")
                else:
                    GTEx_InfoList2[0:len(GTEx_InfoList2)] = ["."]*2

                outsnpline.extend(GTEx_InfoList + GTEx_InfoList2)

                eGeneTSS = ''
                #add eGene-pLI and eGene-TSS
                if egene_GeneName == "--":
                    outsnpline.append("\t".join([".","._."]))
                    eGeneTSS = "."
                else:
                    if dicpLI.get(egene_GeneName,False):
                        outsnpline.append(dicpLI[egene_GeneName])
                    else:
                        outsnpline.append(".")
                        print("eGene genename has no pLI: "+egene_GeneName)
                    if dicTSS.get(egene_GeneName,False):
                        outsnpline.append(dicTSS[egene_GeneName])
                        eGeneTSS = dicTSS[egene_GeneName].split("_")[-1]
                        eGeneTSS_chr = dicTSS[egene_GeneName].split("_")[0]
                    else:
                        outsnpline.append("._.")
                        eGeneTSS = "."
                        print("eGene genename has no TSS: "+egene_GeneName)

                #add the |POS-eGeneTSS|
                length1 = ''
                pos = esnp.split("_")[1]
#                print eGeneTSS
                if not eGeneTSS == ".":
                    if not esnp.split("_")[0] == eGeneTSS_chr:#eSNP chr same as the TSS chr
                        length1 = "not one chr"
                    else:
                        length1 = str(abs(int(pos)-int(eGeneTSS)))
                else:
                    length1 = "."
                outsnpline.append(length1)

                #add eSNP chr,pos,ref,alt etc.
                outsnpline.extend(esnp.split("_")+snp_line[pFDR+1:pgenename+1])

                #add eSNP_genepLI and eSNP_geneTSS
                eSNPGeneTSS = ''
                if esnp_GeneName == "." or "," in esnp_GeneName:
                    outsnpline.append("\t".join(['.','._.']))
                    eSNPGeneTSS = '.'
                else:
                    if dicpLI.get(esnp_GeneName,False):
                        outsnpline.append(dicpLI[esnp_GeneName])
                    else:
                        outsnpline.append(".")
                        print("eSNP genename has no pLI: "+esnp_GeneName)
                    if dicTSS.get(esnp_GeneName,False):
                        outsnpline.append(dicTSS[esnp_GeneName])
                        eSNPGeneTSS = dicTSS[esnp_GeneName].split('_')[-1]
                        eSNPGeneTSS_chr = dicTSS[esnp_GeneName].split('_')[0]
                    else:
                        outsnpline.append("._.")
                        print("eSNP genename has no TSS: "+esnp_GeneName)
                        eSNPGeneTSS = '.'
                length2 = ''
                if not eSNPGeneTSS == ".":
                    if not esnp.split("_")[0] == eSNPGeneTSS_chr:
                        length2 = "not one chr"
                    else:
                        length2 = str(abs(int(pos)-int(eSNPGeneTSS)))
                else:
                    length2 = "."
                outsnpline.append(length2)

                #add other anno info for eSNP
                outsnpline.extend(snp_line[pgenename+1::])
  	        out.write("\t".join(outsnpline)+"\n")
    print("Add the information of Info done!")


import os,sys
import argparse
import re
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = 'Add information for eQTL result(fastQTL) from the file provided .\nContact: lili@novogene.com', formatter_class = RawTextHelpFormatter)
parser.add_argument('--in', metavar = 'File', required = True, help = 'The file which want to be added information. Notice that file must be tab seperated.')
parser.add_argument('--addfile1',metavar = 'String', required = True, help = 'The information file')
parser.add_argument('--addfile2',metavar = 'String', required = True, help = 'The information file')
parser.add_argument('--TSS',metavar = 'String', required = True, help = 'The information file')
parser.add_argument('--pLI',metavar = 'String', required = True, help = 'The information file')
parser.add_argument('--geneinfo',metavar = 'String', required = True, help = 'The information file')
parser.add_argument('--pwd',help="Path of Project Directory for analysis(defalut: ./)",default=os.path.abspath('./'))
parser.add_argument('--out',metavar = 'String', required = True, help = 'Out put prefix.')
parser.add_argument('-v', '-V', help="Show version number and exit.", action='version', version='version: 0.1')

def main():
    argv = vars(parser.parse_args())

    inputfile = argv['in'].strip()
    out =  argv['out'].strip()+'.xls'
#    addfileL = []
#    if ',' in argv['addfile'].strip():
#        addfileL = argv['addfile'].strip().split(",")
#    else:
#        addfileL = argv['addfile'].strip()
    addfile1 = argv['addfile1'].strip()
    addfile2 = argv['addfile2'].strip()
    tssf = argv['TSS'].strip()
    pLIf = argv['pLI'].strip()
    geneinfof = argv['geneinfo'].strip()

    analydir = argv['pwd'].strip()
    analydir = os.path.abspath(analydir)
    if not os.path.exists(analydir):
        print 'The '+analydir+' not exist'
        exit(1)

    dicInfo1 = {}
    dicInfo1 = getInfo(addfile1)
    dicInfo2 = {}
    dicInfo2 = getInfo(addfile2)

    dicgenename = {}
    dicgenename = geteGeneName(geneinfof)
    dicTSS = {}
    dicTSS = getTSS(tssf)
    dicpLI = {}
    dicpLI = getpLI(pLIf)


    addInfo(inputfile,dicInfo1,dicInfo2,dicgenename,dicpLI,dicTSS,out)



if __name__ == '__main__':
    main()
