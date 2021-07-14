'''
refhashes.py: everything known about every reference we've seen.

    key: refhash (sha256 of sorted list of contigs and lengths)

 common: common name, should contain no spaces
    url: authoritative source
   info: description
    md5: md5 of the file at url (or noted in info)
   sha1: sha1 of the file at url (or noted in info)
  build: major nucleotide sequence release
  style: contig naming style for chromosomes

The following all refer to human chromosome 1. The NCBI couldn't settle on a
single convention, and uses different long accession IDs for GenBank and
RefSeq.

     EBI: 1
    UCSC: chr1
 GenBank: CM000663.2
  RefSeq: NC_000001.11

References in the same "build" use the same nucleotide sequences for
chromosomes (with minor difference such as lower case repeat masked regions)
but some build members may contain a subset of contigs of others (eg. decoy
sequences, alt alignments, patches, etc).

References with the same "style" use the same naming convention for the
chromosomes commonly known as 1-22, X, Y, M.

Calls from different references can only be compared if both build and style
match.

'''
from enum import Enum, auto, unique

__all__ = [
    'refhashes', 'refstyle'
]

@unique
class refstyle(Enum):
    ''' Every known contig naming style '''
    UNKNOWN = auto()
    EBI = auto()
    GenBank = auto()
    RefSeq = auto()
    UCSC = auto()

refhashes = {
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855": {
        "common": "no_contigs_present",
        "url": "",
        "info": "refhash found no contigs in the input",
        "md5": "d41d8cd98f00b204e9800998ecf8427e",
        "sha1": "da39a3ee5e6b4b0d3255bfef95601890afd80709",
        "build": "",
        "style": refstyle.UNKNOWN,
    },
    "1e4ef0c15393ae133ad336a36a376bb62e564a43a65892966002e75713282aec": {
        "common": "hs37d5",
        "url": "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
        "info": "GRCh37 primary assembly (chromosomal plus unlocalized and unplaced contigs), the rCRS mitochondrial sequence (AC:NC_012920), Human herpesvirus 4 type 1 (AC:NC_007605) and the concatenated decoy sequences (hs37d5cs.fa.gz)",
        "md5": "a07c7647c4f2e78977068e9a4a31af15",
        "sha1": "7850ab33cf6039b37eafb402b0f9c7abef5d3222",
        "build": "GRCh37",
        "style": refstyle.EBI,
    },
    "9d1184b1f957da7a499793e838a6509626bf772c8f437c0972c25f30fbab9fd7": {
        "common": "human_g1k_v37",
        "url": "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz",
        "info": "1000 Genomes build of GRCh37 from October 2009 including newer version of the MT (NC_012920).",
        "md5": "45f81df94f0408d082363e34a081ed81",
        "sha1": "c59d016293fc41c2f4f08e3f174f6ae297fc307e",
        "build": "GRCh37",
        "style": refstyle.EBI,
    },
    "11cc1f8edf231691097cb007d1cdf603fc53d6019ae2899b6cae23be1ab3393b": {
        "common": "GRCh38_full_analysis_set_plus_decoy_hla",
        "url": "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "info": "Based on GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set, combines HLA sequences with the GRCh38 full analysis set and decoy sequences.",
        "md5": "64b32de2fc934679c16e83a2bc072064",
        "sha1": "efaaea68910ee444b2756062b2ae2b990d5cdb71",
        "build": "GRCh38",
        "style": refstyle.UCSC,
    },
    "1f5faf40c2b1b8715e9df75375cb392117a9c5734fca790e6399d7a50e90ebdd": {
        "common": "grch38",
        "url": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
        "info": "GCA_000001405.15_GRCh38_no_alt_analysis_set released February 2016, combines the human decoy sequences from hs38d1 (GCA_000786075.2) with the GRCh38 'no alt' analysis set.",
        "md5": "a08035b6a6e31780e96a34008ff21bd6",
        "sha1": "70fb7af4dff26bffdf27dbef80caf1f0359d488f",
        "build": "GRCh38",
        "style": refstyle.UCSC,
    },
    "2a1c2512c0d779999183f59350ab7fdb2f89326737af53564d7bd6f1e18aa9be": {
        "common": "hg19",
        "url": "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
        "info": "The UCSC build of the February 2009 GRCh37 assembly.",
        "md5": "806c02398f5ac5da8ffd6da2d1d5d1a9",
        "sha1": "4910f6eea4fd5008b9eded34bddd9a5c9641cb93",
        "build": "GRCh37",
        "style": refstyle.UCSC,
    },
    "37a6fd3d23b48379f4cf17f2fb56faf4cd73226225bc8d7b0f9204e93e3cc901": {
        "common": "hg38",
        "url": "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "info": "The UCSC build of the December 2013 GRCh38 assembly.",
        "md5": "1c9dcaddfa41027f17cd8f7a82c7293b",
        "sha1": "8e8ae3f73d61c3ec8c2477334199557128946276",
        "build": "GRCh38",
        "style": refstyle.UCSC,
    },
    "870d360f9f3ac29b34cfee2a0ab74a78775fcb117ec6c305b243fbf042f7c842": {
        "common": "GRCh38_GenBank",
        "url": "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/",
        "info": "NCBI GRCh38 (GenBank), md5 and sha1 are for GCA_000001405.15_GRCh38_genomic.fna.gz",
        "md5": "9313516997244892bbeceef731fc607e",
        "sha1": "af9f13e9522af59777ca87f2f49891f3de7694f9",
        "build": "GRCh38",
        "style": refstyle.GenBank,
    },
    "65c1965612d479ea8b0530e87a03fe4f9f14dbc11b83274ecdaa700f08d1a1b9": {
        "common": "GRCh38_RefSeq",
        "url": "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/",
        "info": "NCBI GRCh38 (RefSeq), md5 and sha1 are for GCF_000001405.26_GRCh38_genomic.fna.gz",
        "md5": "3d03cc56a6a65413b8abc5a3239d78cd",
        "sha1": "a6d52e51cf7dc5c383bd4d33e8f2970c5afb1e98",
        "build": "GRCh38",
        "style": refstyle.RefSeq,
    },
    "0238e98a614a374caee7b393a2b92f81c714f4075de9c1dce7b9194e044a427c": {
        "common": "GRCh38.p13_GenBank",
        "url": "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/",
        "info": "NCBI GRCh38.p13 (GenBank), md5 and sha1 are for GCA_000001405.28_GRCh38.p13_genomic.fna.gz",
        "md5": "f28b7146e0f30efa58447eceb32620a3",
        "sha1": "6455d4ffe5661564f00521539d6ebb9f628955c4",
        "build": "GRCh38",
        "style": refstyle.GenBank,
    },
    "8450f353b5b093d0520c73e2b1118d5655154a7e417e7c94b5fb4c1c3f6b5110": {
        "common": "GRCh38.p13_RefSeq",
        "url": "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/",
        "info": "NCBI GRCh38.p13 (RefSeq), md5 and sha1 are for GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
        "md5": "84d56a8f8cd75fdde8f60c4e022f9ab7",
        "sha1": "19025e1902ff6c3657e9c846bc141ed323d2a199",
        "build": "GRCh38",
        "style": refstyle.RefSeq,
    },
    "c9f4018aea11b732113e2b2a57f79e6253c7dd9438bc79d39f35dda93c446948": {
        "common": "Homo_sapiens.GRCh38.dna.primary_assembly",
        "url": "http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        "info": "GCA_000001405.28: Ensemble build of GRCh38, including all toplevel sequence regions but excluding haplotypes and patch.",
        "md5": "f58a039f83d8944ceb79cdd16cbda583",
        "sha1": "1c5aab59b9d971cffcc028e4055d985c648b8488",
        "build": "GRCh38",
        "style": refstyle.EBI,
    },
    "7d36a47714a3b3c19d758cafc4c20bc8d319a0bdd06d9d6ab786897c64f2a8a0": {
        "common": "Homo_sapiens.GRCh38.dna.toplevel",
        "url": "http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
        "info": "GCA_000001405.28: Ensemble build of GRCh38, including chromsomes, regions not assembled into chromosomes and N padded haplotype/patch regions.",
        "md5": "ae85c8481141a9399ec7576ad6b0e259",
        "sha1": "a02f246e6682754b80c53e665afab5119d771de2",
        "build": "GRCh38",
        "style": refstyle.EBI,
    },
    "4677415a4ecde9fceca5280f41bd984d90ccd7cd715489ab92fb43250b3cf432": {
        "common": "e_coli_k12_ASM584v1",
        "url": "",
        "info": "Spiral e.coli test reference, md5 and sha1 are for source.fasta.",
        "md5": "d57dd10fbf960a10a05c242094095501",
        "sha1": "edb5e9c75764cbbde097fbbc6bd28408c53cc22a",
        "build": "e_coli",
        "style": refstyle.UNKNOWN,
    },
    "30fad211ba7f62ded67e44167d223332fedb8d8ede00e9bebc8cd633d70664ef": {
        "common": "lambdaToyData",
        "url": "",
        "info": "Spiral lambdaToyData test reference v6.0.0",
        "md5": "9d236d78cfeae6c1b13702ffa3da974f",
        "sha1": "02d0ef55ca950e83dc8efc352ab04a1886037040",
        "build": "lambdaToyData",
        "style": refstyle.UNKNOWN,
    },
    "3d613562b94854287581e37855661ac50e29d109fcbdb73178eca0e91d62012a": {
        "common": "chm13",
        "url": "https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz",
        "info": "Complete T2T reconstruction of a human genome from long reads.",
        "md5": "6d827b6512562630137008830c46e1ac",
        "sha1": "53659e2cec185072f4bf6bb4fb022ec1147c1131",
        "build": "T2T",
        "style": refstyle.UCSC,
    },
    "7847c0ab77e76b50834dbc15bb5d104cacce3bbaa9ec443ae49ab25297312768": {
        "common": "GRCh37.p13_RefSeq",
        "url": "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz",
        "info": "NCBI GRCh37.p13 (RefSeq)",
        "md5": "46e212080d30b1a24abec3eab36dbacd",
        "sha1": "375a4c03856c2989080974a93aba641f09b8cc11",
        "build": "GRCh37",
        "style": refstyle.RefSeq,
    },
}
