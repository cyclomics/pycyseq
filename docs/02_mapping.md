# Mapping

A CySeq read can have multiple mappings to the reference genome. There can be two reasons for this behaviour:

1. The DNA sequence in the CySeq read is not unique in the reference genome, therefore it is not possible to know the origin of the DNA sequence. This problem is not unique to CySeq, but to any sequencing technology.
2. CySeq reads are circular, and it is not possible to know the starting point of a circle. We determine the most likely starting point based on alignment. This approach works for the majority of the reads, but sometimes, there is more than one equally likely starting position.

Below is the description on how these ambiguous cases are annotated for both single-strand and double-strand DNA molecules.

## Single-strand DNA

Single-strand DNA (ssDNA) alignments will be flagged as primary, and can be flagged as supplementary and/or secondary. All the CySeq reads will have a primary alignment, and may have additional supplementary and/or secondary alignments if multiple equally likely alignments exist.

### Unique mapping

If the DNA sequence uniquely maps to a single position in the reference genome it will only have one alignment. This alignment will be flagged as primary.

!!! example

    ```
          Reference:   ...GATCGATGATAGGGTCGATAGCTAGCTGAACGCGATGCAGT...
    Aln 1 (primary):             GATAGGGTCGATAG
    ```

### Unique mapping but multiple rotations

Because CySeq reads originate from circular DNA, sometimes there might me more than one equally likely mapping at the same genomic position. Since we are rotating the sequence of the CySeq read, we call these multiple mappings rotations.

For example, the first `G` in the CySeq read matches the `*`-marked `G` in the reference genome. However, because the CySeq sequence comes from a circle, we can rotate the sequence and therefore the same `G` could also map to the `&`-marked `G` in the reference genome.

The primary alignment will always be the one that aligns to the 5' most position of the reference genome. Supplementary alignments will always rotate bases from the 5' to the 3' relative to the reference genome.

!!! example

    ```
    CySeq sequence: GATAGGGTCGATAG

                                       *             &            
          Reference      :   ...GATCGATGATAGGGTCGATAGGTAGCTGAACGCGATGCAGT...
    Aln 1 (primary)      :             GATAGGGTCGATAG
    Aln 2 (supplementary):              ATAGGGTCGATAGG
    ```

### Non-unique mapping

It is possible that the CySeq DNA sequence is not unique in the reference genome. We can therefore have multiple mappings from the same sequence. 

Since the mappings are at different positions, and these are not due to rotation, additional mappings are flagged as secondary.

!!! example

    ```
    CySeq sequence: GATAGGGTCGATAG

                                                                
          Reference      :   ...GATCGATGATAGGGTCGATAGATAGCTGGATAGGGTCGATAGA...
    Aln 1 (primary)      :             GATAGGGTCGATAG
    Aln 2 (secondary)    :                                  GATAGGGTCGATAG
    ```

It is also possible that there is a secondary mapping by rotating the DNA sequence. Note that in the secondary alignment, the DNA sequence is a rotation of the initial CySeq sequence. Because it is a mapping at a different location, it is also flagged as secondary.

!!! example

    ```
    CySeq sequence: GATAGGGTCGATAG

                                                                
          Reference      :   ...GATCGATGATAGGGTCGATAGATAGCTGTCGATAGAGATAGGG...
    Aln 1 (primary)      :             GATAGGGTCGATAG
    Aln 2 (secondary)    :                                  TCGATAGGATAGGG
    ```

!!! warning "Mapping qualities"
    When a CySeq read has one or more secondary alignments, the mapping quality of ALL the alignments of that read are set to zero, since any alignment could be the correct one. 
    When a CySeq read only has primary alignments and supplementary alignments, the mapping quality is kept as is calculated by the aligner, in this case, minimap2.


## Double-strand DNA

Double-strand DNA (dsDNA) molecules will produce at least two alignments, one for the forward strand, and one for the reverse strand. The two alignments are paired using the R1/R2 flags. The strand that is first detected in the CySeq concatemer will be annotated as R1, the other strand will be annotated as R2.

Similarly to ssDNA, dsDNA alignments will be flagged as primary, and may also be flagged as supplementary if rotations are possible.

## Unique mapping

If the DNA sequences (fwd and rev) uniquely map to a single position in the reference genome, each strand will have its own unique alignment.

!!! example

    ```
              Reference:   ...GATCGCTGATAGGGTCGATAGCTAGCTGAACGCGATGCAGT...
    Aln Fwd 1 (primary):             GATAGGGTCGATAG
    Aln Rev 1 (primary):           GACTATCCCAGCTATCGA   
    ```

## Unique mapping but multiple rotations

Because CySeq reads originate from circular DNA, sometimes there might me more than one equally likely mapping at the same genomic position. In this case, it might not be possible to determine the precise start/end of the two strands because they were connected linearly. In the circle, the 3' end of forward strand is ligated to the 5' start of the reverse strand, and viceversa. If such case arises, we use supplementary alignments to annotate all the possible rotations.

Notice in the example below how the `A` base marked with an `*` on the forward strand could also belong to the reverse strand, since it would be complementary to the `T` base marked with a `&` on the forward strand. It is therefore possible to rotate the `A` base from the 5' end of the forward strand to the 3' end of the reverse strand (Supplementary alignment 2). 

The same logic can be applied to the `T` base marked with an `&` on the forward strand. After rotating the previous `A` base, we can also rotate it from the from the 5' end of the forward strand to the 3' end of the reverse strand (Supplementary alignment 3).

!!! example

    ```                               
                                         *& 
              Reference:      ...GATCGCATATAGGGTCGATAGCTAGCTGAACGCGATGCAGT...
    Aln Fwd 1 (primary):                 ATAGGGTCGATAG
    Aln Rev 1 (primary):                   TCCCAGCTATCGA   

              Reference:      ...GATCGCATATAGGGTCGATAGCTAGCTGAACGCGATGCAGT...
    Aln Fwd 2 (supplementary):            TAGGGTCGATAG
    Aln Rev 2 (supplementary):            ATCCCAGCTATCGA   

              Reference:      ...GATCGCATATAGGGTCGATAGCTAGCTGAACGCGATGCAGT...
    Aln Fwd 3 (supplementary):             AGGGTCGATAG
    Aln Rev 3 (supplementary):           TATCCCAGCTATCGA   
    ```

!!! warning "Merging alignments"
    It is important when analyzing this data to not mix different pairs of alignments. For example, using Aln Fwd 1 with Aln Rev 2 would be incorrect, since the `A` base would be used by the two alignments simultaneously. Alignment pairs are numbered in the `YO` tag. If you read your data using this library, alignment pairs are automatically generated for you. However, if you prefer to read the data yourself, please keep this into consideration.