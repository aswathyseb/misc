# This script reads a bam file and calculate the coverage of each query aganist
# the mapped reference.


import pysam, sys


def parse(bam, chroms):
    bamfile = pysam.AlignmentFile(bam, "rb")

    aln_len, ref_len, covs = dict(), dict(), dict()


    # Reading bam file chromosome by chromosome.
    for chr in open(chroms, "r"):
        chr = chr.strip()

        # Store the length of chromosome
        ref_len[chr] = bamfile.get_reference_length(chr)

        for read in bamfile.fetch(chr):
            if read.is_unmapped: continue

            # Store alignment length for each query:ref pair
            key = ":".join([read.query_name, read.reference_name])
            aln_len[key] = aln_len[key] + read.query_alignment_length if key in aln_len else read.query_alignment_length

    # Calculate coverage
    for k, v in aln_len.items():
        query, ref = k.split(":")
        rlen = ref_len[ref]
        cov = round((v / rlen) * 100, 2)
        covs[k] = cov

    # print
    for k in sorted(aln_len, key=aln_len.get, reverse=True):
        query, ref = k.split(":")
        mapped_len = str(aln_len[k])
        mapped_cov = str(covs[k])
        chrom_len = str(ref_len[ref])
        # if cov < 2.0: continue
        res = "\t".join([query, ref, mapped_len, chrom_len, mapped_cov])
        print(res)


if __name__ == "__main__":
    bam = "test.bam"
    chroms = "test.txt"
    #bam = sys.argv[1]
    #chroms = sys.argv[2]
    parse(bam, chroms)
