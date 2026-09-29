nextflow.enable.types = true

record SampleMeta {
    id: String
    single_end: Boolean
    tool_args: Map<String,String>
}

record Sample {
    id: String
    single_end: Boolean
    reads: List<Path>
    tool_args: Map<String,String>
}

record Alignment {
    id: String
    bam: Path
}

record AlignedSample {
    id: String
    bam: Path
    bai: Path
}
