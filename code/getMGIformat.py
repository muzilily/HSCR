#!/usr/bin/python

import sys,os,re

infile = sys.argv[1]
outfile  = sys.argv[2]

#for R GgenVisR to plot landscape of SNPs
#3 types file format
#MAF: must contain "Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification"
#MGI: must contain "sample, gene_name, trv_type"
#Custom: must contain "sample, gene, variant_class"


def getPos(name,title,name2=''):
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

def getHet(sample):
    if re.match(r'\d+',sample):
        gt = re.search('(\d+)[\/|\|](\d+)',sample)
        if gt.group(1) == gt.group(2) and gt.group(1) != '0':
            return "hom"
        elif gt.group(1) == gt.group(2) and gt.group(1) == "0":
            return "ref"
        elif gt.group(1) != gt.group(2):
            return "het"
        else:
            return "chet"
    elif sample in [".","./."]:
        return 'noGT'

def getMGI(infile, outfile):
    for i in infile:
        if i.startswith("Priority"):
            title = i.strip().split("\t")
            formatp = title.index("FORMAT")
            orirefp = title.index("Ori_REF")
            genenamep = title.index("GeneName")
            funcp = title.index("Func")#exonic or intronic,splicing....
            exonicfuncp = title.index("ExonicFunc") #missense,stopgain...
            aachangep = title.index("AAChange")
            hgmdp = title.index("HGMD_Disease_ID")
            outfile.write("\t".join(["sample","gene_name","trv_type","gt"]+title[1:formatp]+title[hgmdp:hgmdp+2])+"\n")

        else:
            ii = i.strip().split("\t")
            genename = ii[genenamep]

            #valid values for MGI are: nonsense, frame_shift_del, frame_shift_ins, splice_site_del, splice_site_ins, splice_site, nonstop, in_frame_del, in_frame_ins, missense, splice_region_del, splice_region_ins, splice_region, 5_prime_flanking_region, 3_prime_flanking_region, 3_prime_untranslated_region, 5_prime_untranslated_region, rna, intronic, silent, NA
            mut_typeD = {"frameshift deletion":"frame_shift_del", "frameshift insertion":"frame_shift_ins", "stopgain":"nonsense", "splicing":"splice_site", "missense SNV":"missense", "UTR3":"3_prime_flanking_region", "UTR5":"5_prime_flanking_region", "intronic":"intronic", "intergenic":"NA", "synonymous SNV":"silent", "unknown":"NA"}
            trv_type = ''
            func = ii[funcp]
            exonicfunc = ii[exonicfuncp]
            if func in ["exonic","exonic;splicing"]:
                trv_type = mut_typeD[exonicfunc]
            elif func in ["ncRNA_exnoic","ncRNA_exonic;splicing","ncRNA_intronic","ncRNA_splicing"]:
                trv_type = 'rna'
            else:
                if func in mut_typeD.keys():
                    trv_type = mut_typeD[func]
                else:
                    print "the mutype can not be judged! "+func+"_"+exonicfunc
                    trv_type = "NA"

            sample = ''
            for s in range(formatp+1, orirefp):
                s_gt = ii[s]
                if not getHet(s_gt) in ['ref','noGT']:
                #if not getHet(s_gt) in ['']:
                    sample = title[s]
                    outfile.write("\t".join([sample, genename, trv_type, s_gt]+ ii[1:formatp]+ ii[hgmdp:hgmdp+2])+'\n')

with open(infile,'r') as inf, open(outfile,'w') as outf:
    getMGI(inf,outf)
