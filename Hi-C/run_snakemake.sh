#!/bin/bash

if [[ "$@" == "" ]]; then
  set -- "all"
fi

snakemake -j 933 --cluster-config cluster.json --cluster "{cluster.sbatch} -p {cluster.partition} --cpus-per-task {cluster.n} -t {cluster.time} --mem {cluster.mem} {cluster.moreoptions}" --keep-going --latency-wait 600 --restart-times 2 "$@"
