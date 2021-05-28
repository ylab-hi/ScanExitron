VCF files here


```python
import glob

for i in glob.iglob('*.vcf'):
    name = os.path.basename(i).split('.')[0]
    print(f'ScanNeo.py anno -i {i} -o {name}.vep.vcf')
```    
