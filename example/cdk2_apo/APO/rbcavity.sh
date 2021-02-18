#!/bin/sh
#PJM --rsc-list "node=1"
#PJM --rsc-list "elapse=00:10:00"
#PJM --rsc-list "node-mem=12Gi"
#PJM --mpi "shape=1"
#PJM --mpi "proc=1"
#PJM --mpi "rank-map-bychip"
#PJM --mpi assign-online-node
#PJM --mpi "use-rankdir"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ../../../bin/rbcavity %r:./"
#PJM --stgin-dir "rank=* ../../../data/ %r:./data/"
#PJM --stgin-dir "rank=* ../../../data/filters/ %r:./data/filters/"
#PJM --stgin-dir "rank=* ../../../data/pmf/ %r:./data/pmf/"
#PJM --stgin-dir "rank=* ../../../data/pmf/smoothed/ %r:./data/pmf/smoothed/"
#PJM --stgin-dir "rank=* ../../../data/scripts/ %r:./data/scripts/"
#PJM --stgin-dir "rank=* ../../../data/sf/ %r:./data/sf/"
#PJM --stgin "rank=* ./cdk2_apo.prm %r:./"
#PJM --stgin "rank=* ./CDK2_APO.mol2 %r:./"
#PJM --stgout "rank=0 0:./cdk2_apo.as ./"
#PJM -s

. /work/system/Env_base

export RBT_ROOT=.

./rbcavity -r cdk2_apo.prm -was

