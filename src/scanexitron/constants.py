CHRMS: frozenset[str] = frozenset({
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
    "chrX", "chrY", "chrM",
})

NON_MITO_CHRMS: frozenset[str] = CHRMS - {"chrM"}

CHRMS_DICT: dict[str, str] = {
    "1": "chr1", "2": "chr2", "3": "chr3", "4": "chr4", "5": "chr5",
    "6": "chr6", "7": "chr7", "8": "chr8", "9": "chr9", "10": "chr10",
    "11": "chr11", "12": "chr12", "13": "chr13", "14": "chr14", "15": "chr15",
    "16": "chr16", "17": "chr17", "18": "chr18", "19": "chr19", "20": "chr20",
    "21": "chr21", "22": "chr22", "X": "chrX", "Y": "chrY", "MT": "chrM",
}

REVERSE_CHRMS_DICT: dict[str, str] = {v: k for k, v in CHRMS_DICT.items()}

CANONICAL_SPLICE_SITES: frozenset[str] = frozenset({"GT-AG", "GC-AG", "AT-AC"})
