"""reated on Wed Jun 22 17:22:24 2016

@author: vargas
"""
##################################################################################
###################   STRONG SCALING SUBMISSION SCRIPT     #######################
import os

#nlist = [8,12,16,24,32,64,96,102,128,156,196,212]  # Total of 12 jobs
#nnodes=[1,2,1,2,2,4,6,7,8,10,13,14]

for i in range(1,15):
    with open('hadressIG.x','w') as ssf:
            ssf.write('# @ shell=/bin/bash\n')
            ssf.write('#\n')
            ssf.write('# @ error =jobs/jobMRes'+str(i*16)+'.err.$(jobid)\n')
            ssf.write('# @ output =jobs/jobMRes'+str(i*16)+'.out.$(jobid)\n')
            ssf.write('# @ job_type = parallel\n')
            ssf.write('# @ environment = COPY_ALL\n')
            ssf.write('# @ node_usage= shared\n')
            ssf.write('# @ node = '+str(i)+'\n')
            ssf.write('# @ tasks_per_node = 16\n')
            ssf.write('# @ resources =ConsumableCpus(1)\n')
            ssf.write('# @ network.MPI =sn_all,shared,us\n')
            ssf.write('# @ wall_clock_limit =24:00:00\n')
            ssf.write('# @ notification = always\n')
            ssf.write('# @ queue\n')
            ssf.write('\n')
            ssf.write('module load intel/14.0\n')            
            ssf.write('module load python27/python/2.7\n')
            ssf.write('module load python27/scipy/2015.10\n')
            ssf.write('\n')
            #ssf.write('export\n')           
            #ssf.write('PYTHONPATH=${HOME}/espressopp:${HOME}/espressopp/contrib:${PYTHONPATH}\n')
            #ssf.write('export\n')           
            #ssf.write('LD_LIBRARY_PATH=/u/system/SLES11/soft/intel1/14.0/compiler/lib/intel64:${LD_LIBRARY_PATH}\n')
            #ssf.write('\n')
            ssf.write('source /ptmp/gvargas/code/e++MResV1/ESPRC\n')  # other  /ptmp/gvargas/code/espresso2k16Q2GH/espressopp/ESPRC
            ssf.write('\n')
            ssf.write('mpiexec -n '+str(i*16)+' python waterH.py ')
            ssf.close();
    os.system('llsubmit hadressIG.x;')
    os.system('rm hadressIG.x;')
