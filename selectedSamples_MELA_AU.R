dat <- read.table(file = './ICGC_MELA_AU/sample.MELA-AU.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)

uniqDonors <- unique(dat$submitted_donor_id)

selectedSamples <- NULL

for (u in uniqDonors) {
    # Select all samples with that id, then select the first instance 
    selectedSamples <- dat[which(dat$submitted_donor_id == u),][1,]$icgc_sample_id
    selectedSamples <- c(selectedSamples,selectdSamples)
}

write.table(selectedSamples, "./ICGC_MELA_AU/selectedSamples.MELA-AU.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

