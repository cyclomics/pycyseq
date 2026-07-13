# Downstream analysis

Here we provide with some general guidelines into making a CySeq alignment file compatible with traditional genomics downstream analysis.

## Copy Number Analysis

We would recomend to filter the bam file to only contain aligned reads, primary alignments (to avoid multiple counts due to rotations, see the [Mapping section](02_mapping.md) for further information on this topic), only alignments with a mapping quality of >= 20 (you might be more restrictive), discard duplicate reads, and only R1 reads to avoid counting both forward and reverse strands from dsDNA fragments.

```bash title="Example filter CySeq bam for CNV analysis"
# -q 20 only keep mapping quality of 20 or above
# -f 64 only keep R1 reads
# -F 3844 filter out
#       4 unmapped reads
#     256 secondary alignment reads
#     512 failed qc reads
#    1024 duplicated reads
#    2048 supplementary alignment reads
#    3844 = sum

samtools view -b -q 20 -f 64 -F 3844 -o filtered.bam input.bam
```

## Variant calling

We would recomend to filter the bam file to only contain aligned reads, primary alignments (to avoid multiple counts due to rotations, see the [Mapping section](02_mapping.md) for further information on this topic), only alignments with a mapping quality of >= 20 (you might be more restrictive), and to discard duplicate reads.

In this case, we do not filter out only R1 reads since dsDNA might contain variants only present in one of the strands, particularly in the overhang section which is single-stranded.

```bash title="Example filter CySeq bam for variant calling analysis"
# -q 20 only keep mapping quality of 20 or above
# -F 3844 filter out
#       4 unmapped reads
#     256 secondary alignment reads
#     512 failed qc reads
#    1024 duplicated reads
#    2048 supplementary alignment reads
#    3844 = sum

samtools view -b -q 20 -F 3844 -o filtered.bam input.bam
```
