### Ribosome profiling quality control
* Prepare annotaion files
* Get count at each position for each mapping length
* Get mapping length abundance distribution
* Calculate P site by plotting the coverage around TSS and TSE



### Steps to run the pipeline
1. put all bam files into one folder. Here let's say 01_bam
2. In terminal run:
	** python /path/to/p01_prepare_annotation.py -g /full/path/to/gffFile     
	* Right now it only accepts ncbi gff format. -g needs to be gff or gff3 file
3. In terminal run:
	** python /path/to/p02_bam_count.py -b /path/to/folder/that/has/bamFiles -t thread -m max_length -r read_length
	* Here -b is folder instead of file.
	* max_length is the maximum mapping to length to consider. Usually ribosome covers 30 bps, but due to the experiments done, it usually has a wide range. By setting this, it will consider reads ranging from 20 to the setting values.
	* read_length: read length in raw fastq files.
4. P site calibration.  In terminal run:
	** python /path/to/p03_P_site_cal.py -b /path/to/folder/that/has/bamFiles -g /full/path/to/gffFile -t thread
	* parameters are same as before
5. Based on the P site plot, decide the offset for each mapping length when calculating the count. And put the information into a txt file in the same folder with gff file. The file should have two columns, the first one is mapping length, the second one is offset. And the first line of the file should be the header 'offset	len'.
6. 

#### files starting with f** are modules files that has functions or class to parse the files.
#### files starging with m** are parts of analysis pipeline.

* 03_P_site_cal.py: Ribosome usually protects around 30 nts of mRNA sequence, P-site reflects which AA that ribosome is translating. When getting the position of mapped reads, usually we use 5' end or sometime 3'end of the read to indicate where does the read maps to. In this case, we need to calculate how far is the 5/3'end to the P site, which is called P-site offset. In addition, as said before, ribosome usually protects around 30 nts, but it's not exact 30 bps. So here, we extract read counts at positions near tss and tse for each mapping length separately.

