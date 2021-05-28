import glob

out = open('run_ScanNeo_anno.sh','w')
for i in glob.iglob('*.vcf'):
    name = os.path.basename(i).split('.')[0]
    out.write(f'ScanNeo.py anno -i {i} -o {name}.vep.vcf\n')
out.close()
