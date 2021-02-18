#!/bin/sh
#PJM --rsc-list "node=4"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "node-mem=12Gi"
#PJM --mpi "shape=4"
#PJM --mpi "proc=32"
#PJM --mpi "rank-map-bychip"
#PJM --mpi assign-online-node
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ../../../bin/rbdock %r:./"
#PJM --stgin-dir "rank=* ../../../data/ %r:./data/"
#PJM --stgin-dir "rank=* ../../../data/filters/ %r:./data/filters/"
#PJM --stgin-dir "rank=* ../../../data/pmf/ %r:./data/pmf/"
#PJM --stgin-dir "rank=* ../../../data/pmf/smoothed/ %r:./data/pmf/smoothed/"
#PJM --stgin-dir "rank=* ../../../data/scripts/ %r:./data/scripts/"
#PJM --stgin-dir "rank=* ../../../data/sf/ %r:./data/sf/"
#PJM --stgin "rank=* ./cdk2_apo.prm %r:./"
#PJM --stgin "rank=* ./CDK2_APO.mol2 %r:./"
#PJM --stgin "rank=* ./cdk2_apo.as %r:./"
#PJM --stgin "rank=* ../cdk2a2_32.sdf %r:./"
#PJM --stgout "rank=0 0:./apo_docking_act_out.sd ./"
#PJM --stgout "rank=* %r:./proc%r.out ./"
#PJM -s

. /work/system/Env_base

export RBT_ROOT=.

time -p mpiexec -n 32 ./rbdock -r cdk2_apo.prm -p dock.prm -n 32 -i cdk2a2_32.sdf -o apo_docking_act_out

