#!/bin/bash
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module list

## Edit with your job command
ONE_K="/projects/b1213/resources/1kGP/data_raw"
OUTDIR="simulation-files"

mkdir -p ${OUTDIR}

log_message "**** Generate admixture models ****"
cut -f 1,2 ${ONE_K}/integrated_call_samples_v3.20130502.ALL.panel | \
    sed '1d' | sed -e 's/ /\t/g' > ${OUTDIR}/1k_sampleinfo.tsv

# AFR washington
echo -e "100\tAFR_washington\tCEU\tYRI\tPUR" >> ${OUTDIR}/AFR_washington.dat
echo -e "10\t0\t0.34\t0.65\t0.01" >> ${OUTDIR}/AFR_washington.dat

log_message "**** Job ends ****"
