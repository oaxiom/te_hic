
# Note that the number of peaks is too small to be able to call CTCF as a former.

contactZ -i ../data/testbedpe.chr10.bedpe.gz -p ../data/CTCF.chr10.bed.gz -n CTCF
contactZ --gc -i ../data/testbedpe.chr10.bedpe.gz -p ../data/CTCF.chr10.gc.bed.gz -n CTCF_GC

