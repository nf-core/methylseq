
nextflow.preview.types = true

process WRITE_FILE {
    input:
    record(
        name: String,
        items: List<String>,
        newLine: Boolean?
    )

    output:
    file(name)

    exec:
    def path = task.workDir.resolve(name)
    path.delete()
    items.each { item ->
        path << item
        if( newLine )
            path << '\n'
    }
}
