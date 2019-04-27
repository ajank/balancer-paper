 {
	if($0~/\tCO:/){
		if($0~/\tCO:Z:REF/ && $0~/\tCO:Z:ALT/){
			c["both    "]++;
		} else {
			if($0~/\tCO:Z:REF/){
				c["pure_ref"]++;
			} 
			if($0~/\tCO:Z:ALT/){
				c["pure_alt"]++;
			} 
			if($0~/\tCO:Z:UNCLEAR/){
				c["unclear "]++;
			}
		}
	} else {
		c["no_entry "]++;
	}
} 
END {
	OFS="\t";
	for (x in c) {
		print x, c[x];
	}
}
