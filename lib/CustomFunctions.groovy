@groovy.util.logging.Slf4j
class CustomFunctions {

  static void validateParams(workflow, params) {
    // Check that Nextflow version is up to date enough
    // try / throw / catch works for NF versions < 0.25 when this was implemented
    try {
    if( ! workflow.nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
    } catch (all) {
    log.error "====================================================\n" +
                "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
                "  Pipeline execution will continue, but things may break.\n" +
                "  Please run `nextflow self-update` to update Nextflow.\n" +
                "============================================================"
    }
    // Show a big error message if we're running on the base config and an uppmax cluster
    if( workflow.profile == 'standard'){
        if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
            log.error "====================================================\n" +
                    "  WARNING! You are running with the default 'standard'\n" +
                    "  pipeline config profile, which runs on the head node\n" +
                    "  and assumes all software is on the PATH.\n" +
                    "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                    "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                    "============================================================"
        }
    }
  }

  static void workflowReport(summary, workflow, params, baseDir) {

    // Set up the e-mail variables
    def subject = "[nf-core/methylseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/methylseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    email_fields['summary']['Container'] = workflow.container
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/methylseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/methylseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    def ngimethylseqlogo = new File("$baseDir/assets/methylseq_logo.png").bytes.encodeBase64().toString()
    def scilifelablogo = new File("$baseDir/assets/SciLifeLab_logo.png").bytes.encodeBase64().toString()
    def ngilogo = new File("$baseDir/assets/NGI_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:ngimethylseqlogo/, "data:image/png;base64,$ngimethylseqlogo")
    email_html = email_html.replaceAll(~/cid:scilifelablogo/, "data:image/png;base64,$scilifelablogo")
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/methylseq] Pipeline Complete"

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }
  }


}