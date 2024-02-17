Hmel209001o	5125000	5450000
hhec090001 5575000 8975000
hele090001 6500000 11175000
hpar090001 6425000 9625000

## load modules
module load perl/5.26.1-fasrc01
module load perl-modules/5.26.1-fasrc01

## subset to RepeatRegion
head -n3 hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out > hmelv25.bigScaffs.Heliconius/hmelv25.chr09_RepeatRegion.out
head -n3 hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hhec.bigScaffs.Heliconius/hhec.chr09_RepeatRegion.out
head -n3 hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hele.bigScaffs.Heliconius/hele.chr09_RepeatRegion.out
head -n3 hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hpar.bigScaffs.Heliconius/hpar.chr09_RepeatRegion.out

egrep "Hmel209001o" hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out | awk '{ if ($7 >= 5125000 && $6 <= 5450000)  print }' >> hmelv25.bigScaffs.Heliconius/hmelv25.chr09_RepeatRegion.out
egrep "hhec09" hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 5575000 && $6 <= 8975000)  print }' >> hhec.bigScaffs.Heliconius/hhec.chr09_RepeatRegion.out
egrep "hele09" hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 6500000 && $6 <= 11175000) print }' >> hele.bigScaffs.Heliconius/hele.chr09_RepeatRegion.out
egrep "hpar09" hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 6425000 && $6 <= 9625000)  print }' >> hpar.bigScaffs.Heliconius/hpar.chr09_RepeatRegion.out

# get counts
awk '{if ($16 != "*") print }' hmelv25.bigScaffs.Heliconius/hmelv25.chr09_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hmelv25.Rfam.RR09.txt
awk '{if ($16 != "*") print }' hhec.bigScaffs.Heliconius/hhec.chr09_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hhec.Rfam.RR09.txt
awk '{if ($16 != "*") print }' hele.bigScaffs.Heliconius/hele.chr09_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hele.Rfam.RR09.txt
awk '{if ($16 != "*") print }' hpar.bigScaffs.Heliconius/hpar.chr09_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hpar.Rfam.RR09.txt


# Parser Output to generate landscape of Repeats through time [%div]
perl /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/Parsing-RepeatMasker-Outputs/parseRM.pl \
    -i /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/hmelv25.bigScaffs.Heliconius/hmelv25.chr09_RepeatRegion.out -p -l 50,1

perl /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/Parsing-RepeatMasker-Outputs/parseRM.pl \
    -i /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/hhec.bigScaffs.Heliconius/hhec.chr09_RepeatRegion.out -p -l 50,1

perl /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/Parsing-RepeatMasker-Outputs/parseRM.pl \
    -i /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/hele.bigScaffs.Heliconius/hele.chr09_RepeatRegion.out -p -l 50,1

perl /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/Parsing-RepeatMasker-Outputs/parseRM.pl \
    -i /n/holyscratch01/mallet_lab/fseixas/1.pseudo_references/revisions/TEannotation/hpar.bigScaffs.Heliconius/hpar.chr09_RepeatRegion.out -p -l 50,1


## Repeat Region on Chromosome 8 =====
Hmel208001o	3300000	3475000	
hhec080001 3300000 5250000
hele080001 3550000 5300000
hpar080001 3675000 5525000

## subset to RepeatRegion
head -n3 hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out > hmelv25.bigScaffs.Heliconius/hmelv25.chr08_RepeatRegion.out
head -n3 hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hhec.bigScaffs.Heliconius/hhec.chr08_RepeatRegion.out
head -n3 hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hele.bigScaffs.Heliconius/hele.chr08_RepeatRegion.out
head -n3 hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hpar.bigScaffs.Heliconius/hpar.chr08_RepeatRegion.out

#
egrep "Hmel208001o" hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out | awk '{ if ($7 >= 3300000 && $6 <= 3475000)  print }' >> hmelv25.bigScaffs.Heliconius/hmelv25.chr08_RepeatRegion.out
egrep "hhec08" hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 3300000 && $6 <= 5250000)  print }' >> hhec.bigScaffs.Heliconius/hhec.chr08_RepeatRegion.out
egrep "hele08" hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 3550000 && $6 <= 5300000) print }' >> hele.bigScaffs.Heliconius/hele.chr08_RepeatRegion.out
egrep "hpar08" hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 3675000 && $6 <= 5525000)  print }' >> hpar.bigScaffs.Heliconius/hpar.chr08_RepeatRegion.out

# get counts
awk '{if ($16 != "*") print }' hmelv25.bigScaffs.Heliconius/hmelv25.chr08_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hmelv25.Rfam.RR08.txt
awk '{if ($16 != "*") print }' hhec.bigScaffs.Heliconius/hhec.chr08_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hhec.Rfam.RR08.txt
awk '{if ($16 != "*") print }' hele.bigScaffs.Heliconius/hele.chr08_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hele.Rfam.RR08.txt
awk '{if ($16 != "*") print }' hpar.bigScaffs.Heliconius/hpar.chr08_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hpar.Rfam.RR08.txt




## Repeat Region on Chromosome 4 =====
Hmel204001o	5650000	5875000
hhec040001  8000000  8225000
hele040001  4175000 5150000
hpar040001  4125000 5250000

## subset to RepeatRegion
head -n3 hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out > hmelv25.bigScaffs.Heliconius/hmelv25.chr04_RepeatRegion.out
head -n3 hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hhec.bigScaffs.Heliconius/hhec.chr04_RepeatRegion.out
head -n3 hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hele.bigScaffs.Heliconius/hele.chr04_RepeatRegion.out
head -n3 hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hpar.bigScaffs.Heliconius/hpar.chr04_RepeatRegion.out

#
egrep "Hmel204001o" hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out | awk '{ if ($7 >= 5650000 && $6 <= 5875000)  print }' >> hmelv25.bigScaffs.Heliconius/hmelv25.chr04_RepeatRegion.out
egrep "hhec040001" hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 3900000 && $6 <= 5225000)  print }' >> hhec.bigScaffs.Heliconius/hhec.chr04_RepeatRegion.out
egrep "hele040001" hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 4175000 && $6 <= 5150000) print }' >> hele.bigScaffs.Heliconius/hele.chr04_RepeatRegion.out
egrep "hpar040001" hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 4125000 && $6 <= 5250000)  print }' >> hpar.bigScaffs.Heliconius/hpar.chr04_RepeatRegion.out

# get counts
awk '{if ($16 != "*") print }' hmelv25.bigScaffs.Heliconius/hmelv25.chr04_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hmelv25.Rfam.RR08.txt
awk '{if ($16 != "*") print }' hhec.bigScaffs.Heliconius/hhec.chr04_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hhec.Rfam.RR08.txt
awk '{if ($16 != "*") print }' hele.bigScaffs.Heliconius/hele.chr04_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hele.Rfam.RR08.txt
awk '{if ($16 != "*") print }' hpar.bigScaffs.Heliconius/hpar.chr04_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hpar.Rfam.RR08.txt




## Repeat Region on Chromosome 2 =====
Hmel202001o	4075000	4125000
hhec020001 4525000 4875000
hele020001 4900000 5275000
hpar020001 4800000 5150000

## subset to RepeatRegion
head -n3 hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out > hmelv25.bigScaffs.Heliconius/hmelv25.chr02_RepeatRegion.out
head -n3 hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hhec.bigScaffs.Heliconius/hhec.chr02_RepeatRegion.out
head -n3 hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hele.bigScaffs.Heliconius/hele.chr02_RepeatRegion.out
head -n3 hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out > hpar.bigScaffs.Heliconius/hpar.chr02_RepeatRegion.out

#
egrep "Hmel202001o" hmelv25.bigScaffs.Heliconius/hmelv25.bigScaffs.fasta.out | awk '{ if ($7 >= 4075000 && $6 <= 4125000)  print }' >> hmelv25.bigScaffs.Heliconius/hmelv25.chr02_RepeatRegion.out
egrep "hhec020001" hhec.bigScaffs.Heliconius/hhec-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 4525000 && $6 <= 4875000)  print }' >> hhec.bigScaffs.Heliconius/hhec.chr02_RepeatRegion.out
egrep "hele020001" hele.bigScaffs.Heliconius/hele-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 4900000 && $6 <= 5275000) print }' >> hele.bigScaffs.Heliconius/hele.chr02_RepeatRegion.out
egrep "hpar020001" hpar.bigScaffs.Heliconius/hpar-2-hmelv25.HM.100gap.bigScaffs.fasta.out | awk '{ if ($7 >= 4800000 && $6 <= 5150000)  print }' >> hpar.bigScaffs.Heliconius/hpar.chr02_RepeatRegion.out

# get counts
awk '{if ($16 != "*") print }' hmelv25.bigScaffs.Heliconius/hmelv25.chr02_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hmelv25.Rfam.RR02.txt
awk '{if ($16 != "*") print }' hhec.bigScaffs.Heliconius/hhec.chr02_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hhec.Rfam.RR02.txt
awk '{if ($16 != "*") print }' hele.bigScaffs.Heliconius/hele.chr02_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hele.Rfam.RR02.txt
awk '{if ($16 != "*") print }' hpar.bigScaffs.Heliconius/hpar.chr02_RepeatRegion.out | tail -n +4 | awk '{print $11}' | sort | uniq -c | awk '{printf "%s\t%s\n", $2, $1}' > TEcounts/hpar.Rfam.RR02.txt

