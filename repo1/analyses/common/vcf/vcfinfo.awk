$0!~/^#/ {
	split($8,tmp,";");
	for (i in tmp){
		n=split(tmp[i],tmp2,"=");
		if (n==2) {
			INFO[tmp2[1]] = tmp2[2];
		}
	}
	if (NF>9) {
		n=split($9,format,":"); #FORMAT

		for (smp=10; smp<=NF; smp++) {
			m=split($(smp),tmp,":");
			if (m==n) {
				for(i=1;i<=n;i++){
					SAMPLE[(smp-9),format[i]]=tmp[i]
				}	
			}
		}
	}
}
