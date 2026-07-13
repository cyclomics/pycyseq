# API

The `pycyseq` API is a wrapper around `pysam`. 

SAM/BAM files produced by our [workflows](#compatible-workflows) follow the SAM specifications. But because our data has several CySeq-specific characteristic, the API simplifies common processing and manipulation tasks while reducing the risk of errors that can arise from manually parsing BAM files.

## Quickstart

To read the bame file use the `read_cyseq_alignment_file` iterator, which will return an `AlignmentGroup` object at each iteration.

```python title="Code: read a CySeq file"
from pycyseq import read_cyseq_alignment_file

num_ssDNA = 0
num_dsDNA = 0
for alignment_group in read_cyseq_alignment_file(bam_file):
    
    if alignment_group.is_ssdna:
        num_ssDNA += 1
        continue

    if alignment_group.is_dsdna:
        num_dsDNA += 1
        continue

print(f"ssDNA: {num_ssDNA}")
print(f"dsDNA: {num_dsDNA}")
```

The `AlignmentGroup` object will contain one or more alignment objects, all from the same CySeq read. The alignments can can either be `SingleStrandAlignment` or `DoubleStrandAlignment` depending on the type of original molecule (ssDNA or dsDNA).

To understand why one CySeq read can have more than one alignment, please read our [Mapping](../02_mapping.md) section.
More information is available on the respective API documentations of the two classes: [`SingleStrandAlignment`](./05_singlestrand.md) and [`DoubleStrandAlignment`](./04_doublestrand.md).

## Compatible workflows

To generate and process raw CySeq data please refer to one of our workflows:

- CySeq Amplicon: https://github.com/cyclomics/wf-cyseq-amp
- CySeq Whole-Genome: https://github.com/cyclomics/wf-cyseq-wg