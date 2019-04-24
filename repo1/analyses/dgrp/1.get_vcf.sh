VCFpile=/g/korbel/shared/projects/drosophila_balancer/data/variants/SNVs_pileup/wgs.pileup_new.vtnorm.vcf.gz
VCFfree=/g/korbel/shared/projects/drosophila_balancer/data/variants/SNVs2/wgs.freebayes-k.filter.norm.vcf.gz

mkdir vcf 2> /dev/null

# Get VCFs from freebayes
less $VCFfree | vcfbiallelic | awk '/^#/ || ($1~/chr[23][LR]/ && $10~/0\/1/ && $11~/0\/0/)' | bgzip > vcf/bal-spec.vcf.gz
less $VCFfree | vcfbiallelic | awk '/^#/ || ($1~/chr[23][LR]/ && $10~/0\/1/ && $11~/1\/1/)' | bgzip > vcf/vrg-spec.vcf.gz
less $VCFfree | vcfbiallelic | awk '/^#/ || ($1~/chr[23][LR]/ && $10~/1\/1/ && $11~/1\/1/)' | bgzip > vcf/common.vcf.gz

for vcf in vcf/{bal-spec,vrg-spec,common}.vcf.gz;
do tabix $vcf;
done

