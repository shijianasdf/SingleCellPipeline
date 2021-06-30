import os,sys,glob
import multiprocessing as mp
from multiprocessing import Pool, cpu_count

fs = glob.glob('/home/sonic/Data/scRNA_SCN2A/fq_20181016/P*_*')

cmds = []
for f in fs:
	sid = f.split('/')[-1]
	sid1 = sid

	cmd = ' '.join(['time', 'cellranger', 'count',
					'--id=' + sid1,
					'--transcriptome=/home/sonic/Resources/10X/refdata-gex-mm10-2020-A',
					'--fastqs=' + f,
					'--sample=' + sid1,
					'--force-cells=10000',
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

