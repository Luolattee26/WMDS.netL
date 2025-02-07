# This code is for running CPC2 and CPAT



# How to run CPC2 stabdalone version?
cd $CPC_HOME && python ./bin/CPC2.py -i allLNC_fasta.fa -o output.txt

# How to run CPAT
# The CPAT model is trained on human data, so we need to use the human model
# Human_logitModel&Human_Hexamer is provided in the CPAT package
cpat.py -g allLNC_fasta.fa \
-d Human_logitModel.RData \
-x Human_Hexamer.tsv \
-o cpat
