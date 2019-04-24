snakemake -j 8 \
	--cluster-config Snake.cluster.json \
	--cluster "{cluster.sbatch} -p {cluster.partition} --time {cluster.time} --cpus-per-task {cluster.n} --mem {cluster.mem}" \
	--timestamp \
	--keep-going \
	--latency-wait 60 \
	--restart-times 0

