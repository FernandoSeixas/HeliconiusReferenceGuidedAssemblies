for spp in `cat ~/code/heliconius_seixas/1.pseudo_references/spp.list`; do
    echo $spp; 
    perl ~/code/heliconius_seixas/1.pseudo_references/finalizeAssemblies/finalizeGenomes.pl \
        $spp \
        assembliesMEDUSA/$spp-2-hmelv25.HM.100gap.fasta \
        support/hmelv25.chromSense support/hmelv25.scaff-2-chrom \
        Hmel \
        > $spp-2-hmelv25Glued.finalGenome.fasta;
    perl ~/code/heliconius_seixas/1.pseudo_references/finalizeAssemblies/finalizeGenomes.pl \
        $spp \
        assembliesMEDUSA/$spp-2-heradem.HM.100gap.fasta \
        support/heradem.chromSense support/heradem.scaff-2-chrom \
        Herato \
        > $spp-2-herademGlued.finalGenome.fasta;
done
