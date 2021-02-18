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
#PJM --stgin "rank=* ../../bin/rbdock %r:./"
#PJM --stgin-dir "rank=* ../../data/ %r:./data/"
#PJM --stgin-dir "rank=* ../../data/filters/ %r:./data/filters/"
#PJM --stgin-dir "rank=* ../../data/pmf/ %r:./data/pmf/"
#PJM --stgin-dir "rank=* ../../data/pmf/smoothed/ %r:./data/pmf/smoothed/"
#PJM --stgin-dir "rank=* ../../data/scripts/ %r:./data/scripts/"
#PJM --stgin-dir "rank=* ../../data/sf/ %r:./data/sf/"
#PJM --stgin "rank=* ./1sj0_rdock.prm %r:./"
#PJM --stgin "rank=* ./1sj0_rdock.mol2 %r:./"
#PJM --stgin "rank=* ./1sj0_rdock.as %r:./"
#PJM --stgin "rank=* ./1sj0_ligand.sd %r:./"
#PJM --stgout "rank=0 0:./1sj0_docking_out.sd ./"
#PJM --stgout "rank=* %r:./proc%r.out ./"
#PJM -s

. /work/system/Env_base

export RBT_ROOT=.

mpiexec -n 1 ./rbdock -r 1sj0_rdock.prm -p dock.prm -n 10 -i 1sj0_ligand.sd -o 1sj0_docking_out

