```python
import os
import glob

hla = {}
with open('PRAD.info.with.HLA.txt') as f:
    f.readline()
    for line in f:
        l = line.rstrip().split('\t')
        hla[l[0]] = l[-1]

for i in glob.iglob('../VCFs/*.vep.vcf'):
    name = os.path.basename(i).split('.')[0]
    print(f'ScanNeo.py hla -t 16 -i {i} --alleles {hla[name]} --af PSO -e 9 -p /home/user/IEDB/ -o {name}.tsv')
```