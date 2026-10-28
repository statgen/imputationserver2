#!/usr/bin/env nextflow

/*
========================================================================================
statgen/imputationserver2
========================================================================================
Github : https://github.com/statgen/imputationserver2
Author: TOPMed Imputation Server Team
---------------------------
 */

nextflow.enable.types = true

include {
    RefPanelParams ;
    PhasingParams ;
    EncryptionParams ;
    ServiceParams ;
    ImputationParams ;
    AncestryParams ;
    UserParams ;
    PgsParams ;
    PgsCatalogParams
} from './types/params'

include { VcfMetafile ; ImputedChunk } from './types/workflow'

params {
    project: String
    project_id: String? = null
    project_date: String

    files: String
    output: String

    refpanel_yaml: String? = null
    refpanel: RefPanelParams? = null

    phasing: PhasingParams
    allele_frequency_population: String? = null

    password: String? = null
    encryption: EncryptionParams

    service: ServiceParams
    imputation: ImputationParams
    mode: String // imputation | qc-only
    ancestry: AncestryParams? = null
    send_mail: Boolean
    user: UserParams
    build: String

    chunksize: Integer
    min_samples: Integer
    max_samples: Integer

    merge_results: Boolean

    pgs: PgsParams
    pgscatalog: PgsCatalogParams
}

include { INPUT_VALIDATION } from './workflows/input_validation'
include { QUALITY_CONTROL } from './workflows/quality_control'
include { PHASING } from './workflows/phasing'
include { IMPUTATION } from './workflows/imputation'
include { ENCRYPTION } from './workflows/encryption'
include { ANCESTRY_ESTIMATION } from './workflows/ancestry_estimation'
include { PGS_CALCULATION } from './workflows/pgs_calculation'

def error_exit(message: String) {
    println("::error:: ${message}")
    error(message)
}

def success_exit(message: String) {
    println("::message:: ${message}")
    exit(0, message)
}

workflow {

    main:
    def refpanel: RefPanelParams? = params.refpanel
    if (params.refpanel_yaml) {
        refpanel = (RefPanelUtil.loadFromFile(params.refpanel_yaml) as Object) as RefPanelParams
    }
    if (refpanel == null) {
        error_exit('Either refpanel or refpanel_yaml must be provided.')
    }

    //TODO create json validation file
    def phasing_engine = params.phasing.engine
    if (phasing_engine != 'eagle' && phasing_engine != 'no_phasing') {
        error_exit("For phasing, only options 'eagle' or 'no_phasing' are allowed.")
    }

    def encryption_password = params.password ?: PasswordCreator.createPassword()

    def population = params.allele_frequency_population ?: 'off'

    def files = channel.fromPath(params.files)
    files.ifEmpty {
        error_exit('No vcf.gz input files detected.')
    }

    // Find site files from full pattern and make site file pattern relative
    def sites_pattern = "${refpanel.sites}"
    // def sites = "./${file(refpanel.sites).name}"

    def chromosomes: Channel<String> = channel.fromList(1..22)
        .map { chr -> "${chr}" }
        .mix(channel.of('X', 'MT'))

    def site_files: Channel<Path> = chromosomes.flatMap { chr ->
        def sites_file = file(PatternUtil.parse(sites_pattern, [chr: chr]))
        if (!sites_file.exists()) {
            return [] as List<Path>
        }

        def sites_file_index = file(PatternUtil.parse(sites_pattern + ".tbi", [chr: chr]))
        if (!sites_file_index.exists()) {
            error_exit("Missing tabix index for ${sites_file}")
        }

        return [sites_file, sites_file_index]
    }

    if (params.pgs.enabled && (!params.imputation.enabled || params.mode != 'imputation')) {
        error_exit('PGS calculation requires imputation.enabled = true and mode = imputation')
    }

    println("Welcome to ${params.service.name} (${workflow.manifest.version})")

    if (params.imputation.enabled) {
        validation = INPUT_VALIDATION(
            population,
            files,
            refpanel,
            params.imputation,
            params.phasing,
            params.build,
            params.mode,
            params.chunksize,
            params.min_samples,
            params.max_samples,
        )

        qc_metafiles = QUALITY_CONTROL(
            validation.validated_files,
            validation.validation_report,
            site_files.collect(),
            refpanel,
            params.phasing,
            params.build,
            population,
            params.chunksize,
        )

        // check if QC chunks exist in case QC failed
        qc_metafiles.ifEmpty {
            error("QC step failed")
        }

        if (params.mode == 'imputation') {
            def phased: Channel<VcfMetafile> = qc_metafiles

            if (phasing_engine != 'no_phasing') {
                phased = PHASING(qc_metafiles, refpanel, params.phasing, params.service.threads)
            }

            imputation = IMPUTATION(
                refpanel,
                params.imputation,
                phased,
                params.service.threads,
            )

            if (params.merge_results == true) {
                ENCRYPTION(
                    imputation,
                    encryption_password,
                    params.phasing.engine,
                    refpanel,
                    params.imputation,
                    params.encryption,
                )
            }
        }
    }

    def ancestry_enabled = (params.ancestry != null) && params.ancestry.enabled
    def estimated_ancestry: Value<Path?> = null
    if (ancestry_enabled) {
        estimated_ancestry = ANCESTRY_ESTIMATION(files, params.ancestry, params.build)
    }

    if (params.pgs.enabled) {
        PGS_CALCULATION(imputation, estimated_ancestry, params.pgs, params.pgscatalog)
    }

    onComplete:
    //TODO: use templates
    //TODO: move in EmailHelper class
    if (!workflow.success) {
        def statusMessage = (workflow.exitStatus != null || workflow.errorReport == "QC step failed")
            ? "failed"
            : "canceled"

        def statusText = (statusMessage == "failed")
            ? "Your job failed."
            : "Your job has been canceled."

        if (params.send_mail && params.user.email != null) {
            sendMail(
                to: "${params.user.email}",
                subject: "[${params.service.name}] Job ${params.project} ${statusMessage}",
                body: "Dear ${params.user.name},\n\n${statusText}\n\nMore details can be found at the following link: ${params.service.url}/index.html#!jobs/${params.project_id}",
            )
        }

        error_exit("Imputation job ${statusMessage}.")
    }

    //submit counters for successful imputation jobs
    if (params.mode == 'imputation') {
        println("::submit-counter name=samples::")
        println("::submit-counter name=variants::")
        println("::submit-counter name=input-chunks::")
        println("::submit-counter name=remaining-chunks::")
        println("::submit-counter name=runs::")
        println("::submit-counter name=chromosomes::")
        println("::submit-counter name=chromosomes-per-sample::")
        println("::submit-counter name=liftover::")

        println("::submit-value name=input-phase::")
        println("::submit-value name=chromosome-list::")

        println("::set-value-and-submit name=reference_panel::${refpanel.id}")
        println("::set-value-and-submit name=phasing_engine::${phasing_engine}")
        println("::set-value-and-submit name=genome_build::${params.build}")
    }

    // imputation job
    if (params.merge_results == true && params.encryption.enabled == true) {
        if (params.send_mail && params.user.email != null) {
            sendMail(
                to: "${params.user.email}",
                subject: "[${params.service.name}] Job ${params.project} is complete",
                body: "Dear ${params.user.name},\n\nYour imputation job has finished succesfully. The password for the imputation results is: ${encryption_password}\n\nYou can download the results from the following link: ${params.service.url}/index.html#!jobs/${params.project_id}",
            )
            success_exit("Data have been exported successfully. We have sent a notification email to <b>${params.user.email}</b>")
        }
        else {
            success_exit("Data have been exported successfully. We encrypted the results with the following password: <b>${encryption_password}</b>")
        }
    }
    else {
        //PGS job
        if (params.send_mail && params.user.email != null) {
            sendMail(
                to: "${params.user.email}",
                subject: "[${params.service.name}] Job ${params.project} is complete",
                body: "Dear ${params.user.name},\n\nYour PGS job has finished successfully.\n\nYou can download the results from the following link: ${params.service.url}/index.html#!jobs/${params.project_id}",
            )
            success_exit("Data have been exported successfully. We have sent a notification email to <b>${params.user.email}</b>")
        }
        else {
            success_exit("Data have been exported successfully.")
        }
    }
}
