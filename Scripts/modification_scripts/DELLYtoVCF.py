# This script is used to convert the output of delly (WHOLE GENOME) into the standard vcf format!!
import sys
f=open(sys.argv[1],"r")
lines=f.readlines()
software_result=[]
for x in lines:
  if (x[0].isdigit() or x[0]=="Y" or x[0]=="X"):
                software_result.append(x.split())
f.close()
bmark = [0] * (len(software_result))

with open(sys.argv[2],'w') as f1:
    f1.write('##fileformat=VCFv4.2\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant detected">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">\n#CHROM     POS     ID     REF    ALT     QUAL  FILTER  INFO\n')
    for j in range(len(software_result)):
        if(software_result[j][4]=="<DEL>"):
            f1.write(software_result[j][0]) #CHROM
            f1.write('\t') # tab
            f1.write(software_result[j][1]) #POS
            f1.write('\t'+'.'+'\t'+'.'+'\t'+'<DEL>'+'\t'+'.'+'\t'+'PASS'+'\t'+'SVTYPE=DEL;SVLEN=')
            f1.write(str( int((software_result[j][7].split(";")[4]).split("=")[1] )- int (software_result[j][1] )  ) ) #start pos - end pos
            f1.write(";END=")
            f1.write( (software_result[j][7].split(";")[4]).split("=")[1] ) #END in FORMAT
            f1.write('\n')
