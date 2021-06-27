import os, sys, glob
import multiprocessing as mp
from multiprocessing import Pool, cpu_count

fs = glob.glob("/data1/NGR_SNU_SingleCell/Data/fastq.gz/*")

cmds = []

for f in fs:
    sid = f.split('/')[-1]
    sid1 = sid
    cmd = ' '.join(['time', 'cellranger', 'count', 
    	'--id=' + sid1, 
    	'--transcriptome = /data1/NGR_SNU_SingleCell/Data/refdata-gex-mm10-2020-A',
    	'--fastqs=' + f,
		'--sample=' + sid1,
		'--expect-cells=10000',
		'--localcores=20',
		'--localmem=256'])
    print (cmd)
    cmds.append(cmd)
    # os.system(cmd)

## Creating a pool for parallel processing 
number_threads = 4
pool = mp.Pool(number_threads)
pool.imap_unordered(os.system, cmds) 
pool.close() 
pool.join() 

