# This code is used in our analysis process to convert standard gtf files into exon as well as transcript bed files
# The transcript bed file is used in the next step of conservativeness analysis



# gtf to bed by bash script
cat gencode.v22.annotation.gtf | grep 'exon' | cut -f1,4,5,9 | cut -f1 -d";" |  awk '{print $1, $2-1, $3, $5}' | sed -e 's/ /\t/g' | sed -e 's/\"//g' > all.exon.bed

# gtf to transcripts bed file
cat gencode.v22.annotation.gtf | awk -F '[\t *;]' '/^chr/{if($3=="transcript"){print $1,$4-1,$5,$10,$13,$22,$7,$3}}' OFS="\t" > all.transcripts.bed