#!/bin/bash
mkdir logs 2> /dev/null
snakemake -j 100\
    --cluster-config Snake.cluster.json \
    --cluster "{cluster.sbatch} -p {cluster.partition} --time {cluster.time} --cpus-per-task {cluster.n} --mem {cluster.mem}" \
	--latency-wait 200\
	--timestamp\
	--restart-times 1\
	--keep-going

