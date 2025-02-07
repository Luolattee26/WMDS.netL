# This code was used in the analysis flow of our paper to download lncRNA interactions data from starBase



# miRNA-LNC
curl 'https://rna.sysu.edu.cn/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA=all&clipExpNum=5&degraExpNum=1&pancancerNum=10&programNum=5&program=None&target=all&cellType=all' > starBase_miRNA_lnc.txt
# RNA-RNA(failed)
curl 'https://rna.sysu.edu.cn/encori/api/RNARNA/?assembly=hg38&geneType=lncRNA&RNA=all&interNum=1&expNum=1&cellType=all' > starBase_RNA_lnc.txt
# RBP target
curl 'https://rna.sysu.edu.cn/encori/api/RBPTarget/?assembly=hg38&geneType=lncRNA&RBP=all&clipExpNum=5&pancancerNum=0&target=all&cellType=all' > starBase_RBP_lnc.txt