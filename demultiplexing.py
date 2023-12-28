##gunzip the R1 and R2 AbSeq+SMK File
gunzip B1_Abseq_S1_L001_R*
##Fetch the SMK reads from Abseq_R2 file by fetching 1 line before and 2 lines after the match
##save the following script as test.awk
BEGIN {
    for (i=1; i<=length(seq); i++) {
        regexp = regexp sep substr(seq,1,i-1) "." substr(seq,i+1)
        sep = "|"
    }
}
{ rec = rec $0 ORS }
!(NR % 4) {
    if (rec ~ regexp) {
        printf "%s", rec
    }
    rec = ""
}
awk -v seq="SampleTagX_Sequence" -f test.awk B1STX_SMK_R2.fastq
##For the samplewise SMK R2 file, fetch the corresponding R1 file from the Batch level Abseq_R1 file B1_Abseq_S1_L001_R1_001.fastq
fastq_pair B1_Abseq_S1_L001_R1_001.fastq B1STX_SMK_R2.fastq 
##this will create 4 files, a set of paired fastq files, and a set of un-matched fastq files. Size of  R2 paired.fq file should match with B1STX_SMK_R2.fastq. Remove the single.fq files to save space
rm -rf *.single.fq
mv B1_Abseq_S1_L001_R1_001.fastq.paired.fq B1STX_SMK_R1_paired.fq
mv B1STX_SMK_R2.fastq.paired.fq B1STX_SMK_R2_paired.fq
##confirm the fastq structure is correct, by counting the number of "@" and "+". The numbers should match.
grep "@A" B1STX_SMK_R1_paired.fq | wc -l
grep "+" B1STX_SMK_R1_paired.fq | wc -l
##Use UMI-tools to fetch the unique cell barcodes from the R1 of the SMK reads.
umi_tools whitelist --stdin dir/B1STX_SMK_R2_paired.fq \
					--extract-method=regex \
					--bc-pattern="(?P<cell_1>.{9})(?P<discard_1>GTGA)(?P<cell_2>.{9})(?P<discard_1>GACA)(?P<cell_3>.{9})(?P<umi_1>.{6,9})T{3}.*"  \
					--expect-cells=20000 \
					--knee-method=density \
					--log2stderr > whitelist.txt;
##this may give some error with "no local minima" if the CB diversity is less, use --set-cell-number in that case instead of --expect-cells & also remove --knee-method=density
##this will create a whitelist of cell barcodes (CB) from the SMK R1 file, with an expected cells of 20000, however, if the cell barcodes (CB) count is more/less that that, it'll generate the actual number.
##next, use the cell barcode (CB) whitelist to demultiplex the samplewise fastq file from the combined fastq file.
umi_tools extract	--extract-method=regex \
					--bc-pattern="(?P<cell_1>.{9})(?P<discard_1>GTGA)(?P<cell_2>.{9})(?P<discard_1>GACA)(?P<cell_3>.{9})(?P<umi_1>.{6,9})T{3}.*" \
					--stdin b1_combined_R1.fastq.gz \
					--stdout b1_st1_combined_R1.fastq \
					--read2-in b1_combined_R2.fastq.gz \
					--read2-out=b1_st1_combined_R2.fastq \
					--whitelist=whitelist.txt;
##this will extract the sample wise R1 and R2 file from the combined batch level fastq file, with CB and UMI moved from the read sequence to the read name of both R1 and R2. 
##this will also leave the extracted R1 file without any actual CB and UMI reads, since those have now been moved the fastq header. While only R2 is sufficient for mapping using UMI-tools,
##separate R1 file can be recreated using the following command for use with other mapping tools.
cat b1_st1_combined_R2.fastq | sed -e 's/_.* / /' > b1_st1_combined_edited_R2.fastq ##this will remove the CB and UMI sequence added by the umi tools from the header
##fetch the R1 corresponding to the b1_st1_combined_edited_R2.fastq from the combined fastq file
gunzip -k b1_combined_R1.fastq.gz
fastq_pair b1_combined_R1.fastq.gz b1_st1_combined_edited_R2.fastq
rm -rf *.single.fq
mv b1_combined_R1.fastq.paired.fq b1_st1_R1.fastq
mv b1_st1_combined_edited_R2.fastq.paired.fq b1_st1_R2.fastq
##once done, transfer the files to SevenBridges for analysis, or run the CWL-based BD Rhapsody pipeline on the local server (recommended)
