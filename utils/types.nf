nextflow.enable.types = true

record Sample {
    id: String
    meta: Record
    reads: List<Path>
}

record Alignment {
    id: String
    meta: Record
    bam: Path
}

record AlignedSample {
    id: String
    meta: Record
    bam: Path
    bai: Path
}
