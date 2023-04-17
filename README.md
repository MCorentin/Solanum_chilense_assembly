# Solanum chilense assembly pipeline

Scripts and files used to perform the de novo assembly of Solanum chilense (LA1972).

The gene models predicted by Augustus are available in [S.chilense.LA1972.gff](Augustus/S.chilense.LA1972.gff)

The organelle assemblies and annotations are available under [organelles/](organelles/)

## How to cite

Coming soon

## Table of Contents

- [Solanum chilense assembly pipeline](#solanum-chilense-assembly-pipeline)
- [How to cite](#how-to-cite)
- [Table of Contents](#table-of-contents)
- [Data](#Data)
- [The assembly pipeline](#the-assembly-pipeline)
	- [Versions](#versions)
	- [MaSuRCA](#masurca)
	- [Redundans](#redundans)
	- [SSPACE](#sspace)
	- [Hybrid Scaffold](#hybrid-scaffold)
	- [Arcs + Link](#arcs-+-link)
	- [Pilon](#pilon)
	- [BBMap dedupe](#bbmap-dedupe)
	- [GapFiller](#gapfiller)
	- [ALLHiC](#allhic)
	- [JuiceBox](#juicebox)
- [Running Augustus](#running-augustus)
- [Circos](#circos)
	
## Data

The assembly and reads have been submitted to the Sequence Read Archive under the Bioproject:

 - Bioproject PRJNA880259, link coming soon.

 - The Bionano optical mapping data is available as supplementary file of the Bioproject: "SUPPF_0000004381".


## The assembly pipeline

![Solanum chilense LA1972 assembly pipeline](./figures/chilense_assembly_pipeline_v2.png?raw=true)

### Versions

<details>

<summary>Version of the tools used to produce the assembly</summary>

- MaSuRCA			v3.2.7
- redundans			v0.13a
- SSPACE			v1-1
- Hybrid Scaffold	v1.0
- RefAligner		v1.0
- Arcs				v8.25
- Links				v1.8.6
- Pilon				v1.22
- BBMap				v37.72
- GapFiller			v1-10
- ALLHiC			v0.9.13
- JuiceBox			v1.11.08
- bwa				v0.7.17
- perl				v5.14.4
- samtools			v1.9
- LongRanger		v2.2.2
- Trimmomatic		v0.39
- blast+			v2.6.0
- bedtools			v2.27.1
- LACHESIS 			(git commit: https://github.com/shendurelab/LACHESIS/commit/2e27abb127b037f87982021e86a45f289cc5e3be)
- sra-cleaning		v0.2.0

</details>

### MaSuRCA

Generating the initial assembly from the Illumina and PacBio reads.

````
# First generate the "assemble.sh" script from the config file
masurca -g config.txt

# Then run the assemble.sh script
sh assemble.sh
````

### Redundans

To remove duplications from the contig assembly.

````
python2.7 redundans.py -v -i /path/to/1490_R1.fastq /path/to/1490_R2.fastq /path/to/1494_R1.fastq /path/to/1494_R2.fastq -f final.genome.scf.fasta  -o masurca237_redundans -t 60 --log redundans_masurca237.log
````

### SSPACE

To scaffold the contigs, using information from the long reads.

````
perl SSPACE-LongRead.pl -c dedup.genome.scf.fasta -p Sequelf1000_RSIIALL.fasta -b ./ - t 40
````

### Hybrid Scaffold

Scaffolding using the BioNano optical maps.

````
perl hybridScaffold.pl -B 1 -N 1 -f -n scaffolds.fasta -b exp_refineFinal1_contigs_C.cmap -c hybridScaffold_config.xml -o /path/to/3_hybridScaffold/ -r /path/to/RefAligner
````

Reintegrating the unmapped scaffolds with "HybridScaffold_finish.pl" (see: https://github.com/i5K-KINBRE-script-share/Irys-scaffolding)

````
perl hybridScaffold_finish_fasta.pl -x ../Hybrid_scaffold/hybrid_scaffolds/exp_refineFinal1_contigs_C_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap -s ../Hybrid_scaffold/hybrid_scaffolds/exp_refineFinal1_contigs_C_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta -f ../Hybrid_scaffold/hybrid_scaffolds/scaffolds.fasta
````

### Arcs + Link

Scaffolding with the Chromium 10x data.

Interleaving the chromium 10x reads with LongRanger.
*/path/to/chromium_10x/* is the folder containing the paired 10x fastq files.

````
# First generate the interleaved file
longranger basic --id=chilense_10x_basic --fastqs=/path/to/chromium_10x/

# Then add the barcode to the file (needed for Arcs)
gunzip -c barcoded.fastq.gz | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print"$_\n" if($flag);}' > chromium_interleaved.fastq
````

Aligning the interleaved fastq to the assembly:
````
bwa index scaffolds_genome_post_HYBRID_SCAFFOLD.fasta

bwa mem -t 60 scaffolds_genome_post_HYBRID_SCAFFOLD.fasta -p chromium_interleaved.fastq > chilense_masurca327_reduced_SSPACE_HSf.sam
samtools view -b chilense_masurca327_reduced_SSPACE_HSf.sam -@ 40 > chilense_masurca327_reduced_SSPACE_HSf.bam
# The alignment.fof file will be used as input to Arcs:
echo "/path/to/chilense_masurca327_reduced_SSPACE_HSf.bam" > alignment.fof
````

Launching Arcs:
````
arcs -f scaffolds_genome_post_HYBRID_SCAFFOLD.fasta -a alignments.fof -s 95 -c 5 -l 0 -z 500 -m 30-10000 -d 0 -e 30000 -r 0.05 -v 1
````

Launching Links (*empty.fof* is an empty file).
````
# Creating the checkpoint.tsv file for Links
makeTSVfile.py scaffolds_genome_post_HYBRID_SCAFFOLD.fasta.scaff_s95_c5_l0_d0_e30000_r0.05_original.gv scaffolds_genome_post_HYBRID_SCAFFOLD.Arcs.tigpair_checkpoint.tsv scaffolds_genome_post_HYBRID_SCAFFOLD.fasta

# Launching Links:
LINKS -f ../scaffolds_genome_post_HYBRID_SCAFFOLD.fasta -b ../scaffolds_genome_post_HYBRID_SCAFFOLD.Arcs -s empty.fof -k 20 -l 5 -t 2 -v 1
````

### Pilon

Correcting errors and misassemblies with the Illumina reads (repeated twice).

````
# First, rename the sequence names to avoid errors related to ID length
sed 's/,/_/g' scaffolds_genome_post_HYBRID_SCAFFOLD.Arcs.scaffolds.fa > chilense_masurca327_scaffoldsReduced_sspace_HSfinish_RENAMED.fasta

bwa index chilense_masurca327_scaffoldsReduced_sspace_HSfinish_RENAMED.fasta

bwa mem -t 60 chilense_masurca327_scaffoldsReduced_sspace_HSfinish_RENAMED.fasta /path/to/1490_R1.fastq /path/to/1490_R2.fastq > chilense_1490.sam
bwa mem -t 60 chilense_masurca327_scaffoldsReduced_sspace_HSfinish_RENAMED.fasta /path/to/1494_R1.fastq /path/to/1494_R2.fastq > chilense_1494.sam

samtools view -b -@ 80 chilense_1490.sam > chilense_1490.bam
samtools view -b -@ 80 chilense_1494.sam > chilense_1494.bam

samtools sort -@ 80 chilense_1490.bam > chilense_1490.sorted.bam
samtools sort -@ 80 chilense_1494.bam > chilense_1494.sorted.bam

samtools index -@ 80 chilense_1490.sorted.bam
samtools index -@ 80 chilense_1494.sorted.bam
````

Running pilon:
````
java -jar -Xmx750G pilon-1.22.jar --genome chilense_masurca327_scaffoldsReduced_sspace_HSfinish_RENAMED.fasta --frags chilense_1490.sorted.bam --frags chilense_1494.sorted.bam --output pilon_corrected_assembly.fasta --outdir ./ --changes --fix all --threads 60
````

### BBmap dedupe

Removing duplicated contigs.

````
dedupe.sh in=pilon_corrected_assembly.fasta out=chilense_pilonRound2_deduped.fa outd=duplicateScaffolds.fasta threads=60 storequality=f absorbrc=t touppercase=t minidentity=90 minlengthpercent=0 minoverlappercent=0 maxsubs=40000 maxedits=5000 minoverlap=1000 k=31 -eoom -Xmx300G
````

### GapFiller

Filling the gaps (stretch of N's) with the Illumina paired-end reads.

````
/usr/bin/perl5.22.1 GapFiller.pl -l libraries.txt -s chilense_pilonRound2_deduped.fa -m 20 -T 30 -b chilense_pilonr2_bbmap1
````

### ALLHiC

Hi-C reads were used to first correct the assembly, and then order and orient the scaffolds into chromosomes.

First trim the Hi-C reads:

````
java -jar trimmomatic-0.39.jar PE /path/to/Wild_Tomato_R1.fastq.gz /path/to/Wild_Tomato_R2.fastq.gz \
 Wild_Tomato_R1_paired.fq.gz Wild_Tomato_R1_unpaired.fq.gz \
 Wild_Tomato_R2_paired.fq.gz Wild_Tomato_R2_unpaired.fq.gz \
 SLIDINGWINDOW:4:20 MINLEN:50 -threads 60 -trimlog ./trim.log
````

Align the Hi-C reads to the assembly:

````
bwa index -a bwtsw -p chilense_pilonr2_bbmap1.gapfilled.final.RENAMED chilense_pilonr2_bbmap1.gapfilled.final.RENAMED.fa

# -SP   Align the pairs as independent single-end reads but still with all 
#       pair-related flags added properly.
# -5    for split alignment, take the alignment with the smallest coordinate as primary.
# -F 2316:  only include reads with NONE of the FLAGS in INT (2316):
#           read unmapped (0x4)
#           mate unmapped (0x8)*
#           not primary alignment (0x100)
#           supplementary alignment (0x800)
bwa mem -t 50 -SP -5 chilense_pilonr2_bbmap1.gapfilled.final.RENAMED Wild_Tomato_R1_paired.fq.gz Wild_Tomato_R2_paired.fq.gz | samtools view -h -b -F 2316 > Wild_Tomato_HiC.bam

samtools sort -n Wild_Tomato_HiC.bam -o Wild_Tomato_HiC.sorted.bam -@ 50

# First pre-process the BAM file with the perl script from LACHESIS to remove noise
perl ./LACHESIS/src/bin/PreprocessSAMs.pl Wild_Tomato_HiC.sorted.bam chilense_pilonr2_bbmap1.gapfilled.final.RENAMED.fa
samtools index Wild_Tomato_HiC.sorted.REduced.sorted.bam
````

Correct the assembly based on the alignment:

````
ALLHiC_corrector -m Wild_Tomato_HiC.sorted.REduced.sorted.bam -r chilense_pilonr2_bbmap1.gapfilled.final.RENAMED.fa	-o chilense.corrected.fasta -t 60 > Correct.log
````

Align the Hi-C reads against the corrected assembly:

````
bwa index -a bwtsw -p chilense.corrected chilense.corrected.fasta

bwa mem -t 50 -SP -5 chilense.corrected Wild_Tomato_R1_paired.fq.gz Wild_Tomato_R2_paired.fq.gz | samtools view -h -b -F 2316 > chilense_corrected_HiC.bam
samtools sort -n chilense_corrected_HiC.bam -o chilense_corrected_HiC.sorted.bam -@ 50
# As before, preprocess the BAM file 
perl ./LACHESIS/src/bin/PreprocessSAMs.pl chilense_corrected_HiC.sorted.bam chilense.corrected.fasta
````

Run ALLHiC:

````
# Filter the alignments with quality < 40 (to keep best alignments)
samtools view -b -q 40 -@ 50 chilense_corrected_HiC.sorted.REduced.bam > chilense_corrected_HiC.sorted.REduced.unique.bam

# Partition the scaffolds into clusters:
ALLHiC_partition -b chilense_corrected_HiC.sorted.REduced.unique.bam -r chilense.corrected.fasta -e MBOI -k 12 -m 25

# Extract the CLM files:
allhic extract chilense_corrected_HiC.sorted.REduced.unique.bam chilense.corrected.fasta --RE MBOI

# Optimize the clusters:
for cluster in $(find ./ -name "chilense_corrected_HiC.sorted.REduced.unique.counts_GATC.12g?*.txt");
do
    cmd="allhic optimize ${cluster} chilense_corrected_HiC.sorted.REduced.unique.clm &";
    eval ${cmd}
done

# Build the assembly
ALLHiC_build chilense.corrected.fasta
````

### JuiceBox

Create the needed .hic and .assembly files:

````
git clone https://github.com/phasegenomics/juicebox_scripts.git
git clone --recursive https://github.com/phasegenomics/matlock.git matlock ; cd matlock ; make

cd ../

./juicebox_scripts/juicebox_scripts/agp2assembly.py groups.agp groups.assembly

# BAM should be sorted by read name (-n)
./matlock/bin/matlock bam2 juicer chilense_corrected_HiC.sorted.bam out.links.txt 
sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
bash ./3d-dna/visualize/run-assembly -visualizer.sh groups.assembly out.sorted.links.txt
````

Here, use JuiceBox to correct misjoins from the .hic and .assembly files.

Then recreate the assembly from the reviewed .assembly file:
````
python ./juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py -a groups.reviewed.assembly -f ../2_correct_assembly/chilense.corrected.fasta 
````

## Cleaning after the SRA report

Removing and trimming the sequences listed in the Contamination.txt file received from the Sequence Read Archive after submission:

````
python sra-cleaning.py -a chilense.final.corrected.renamed.filtered.fasta -g augustus.with.hints.filtered.gff -c My_Contamination.txt -o chilense_cleaned/
````

The sra-cleaning.py script is available at https://github.com/MCorentin/sra-cleaning 


### Running Augustus

First, repeats are masked with ReapeatMasker:
````
RepeatMasker -pa 50 --noisy --xsmall --lib repeats_master.fasta chilense.final.corrected.renamed.filtered.fasta
````

The hints are created with the "bam2hints.pl" script from Augustus. "chilense.rna.merged.bam" comes from the
alignment of the RNA-seq to the assembly with STAR. The resulting bam files are merged with "samtools merge".
````
/path/to/Augustus/auxprogs/bam2hints/bam2hints --intronsonly --in=chilense.rna.merged.bam --out=introns_chilense.gff
````

Runing Augustus with the hints:
````
augustus --species=tomato --UTR=on --softmasking=on --extrinsicCfgFile=extrinsic.M.RM.E.W.P.tomato.cfg --hintsfile=introns_chilense.gff --allow_hinted_splicesites=atac --alternatives-from-evidence=on chilense.final.corrected.renamed.filtered.fasta.masked > augustus.hints.gff
````

Extracting the sequences from the gff:
````
/path/to/Augustus/scripts/getAnnoFasta.pl augustus.hints.gff --seqfile chilense.final.corrected.renamed.filtered.fasta.masked
````

Output:
- `augustus.hints.aa`           =  Amino acid sequences
- `augustus.hints.cdsexons`     =  Coding exon positions on genome
- `augustus.hints.codingseq`    =  Coding sequences only
- `augustus.hints.mrna`         =  mRNA sequences

### Circos

Below is a circos plot representing the final assembly.

From the outside to the inside, each layer represents:
- The list of pseudomolecules
- Gene density (purple = low density, yellow = high density)
- SNP density against S. lycopersicum
- SNP density against S. pennellii
- SNP density against S. chilense LA3111
- GC content (red = lower than genome mean, green = higher)

Please refer to our publication for more details.

![Solanum chilense LA1972 circos](./figures/circos.png?raw=true)
