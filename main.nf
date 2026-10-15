#!/usr/bin/env nextflow

/*
========================================================================================
statgen/imputationserver2
========================================================================================
Github : https://github.com/statgen/imputationserver2
Author: TOPMed Imputation Server Team
---------------------------
 */

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

workflow {

    main:
    refpanel = params.refpanel_yaml ? (record(RefPanelUtil.loadFromFile(params.refpanel_yaml)) as RefPanelParams) : params.refpanel
    if (refpanel == null) {
        exit(1, 'Either refpanel or refpanel_yaml must be provided.')
    }

    //TODO create json validation file
    def phasing_engine = params.phasing.engine

    if (phasing_engine != 'eagle' && phasing_engine != 'no_phasing') {
        println("::error:: For phasing, only options 'eagle' or 'no_phasing' are allowed.")
        exit(1)
    }

    encryption_password = params.password ?: PasswordCreator.createPassword()

    population = params.allele_frequency_population ?: 'off'

    files = channel.fromPath(params.files)
    files.ifEmpty {
        println('::error:: No vcf.gz input files detected.')
        exit(1)
    }

    // Find site files from full pattern and make site file pattern relative
    sites_pattern = "${refpanel.sites}"
    sites = "./${file(refpanel.sites).fileName}"

    site_files_ch = channel.of(1..22, 'X', 'MT')
        .map { it ->
            def sites_file = file(PatternUtil.parse(sites_pattern, [chr: it]))
            def sites_file_index = file(PatternUtil.parse(sites_pattern + ".tbi", [chr: it]))

            if (!sites_file.exists()) {
                return null
            }

            if (sites_file.exists() && !sites_file_index.exists()) {
                error("Missing tabix index for " + sites_file)
            }

            return tuple(sites_file, sites_file_index)
        }

    println("Welcome to ${params.service.name} (${workflow.manifest.version})")

    if (params.imputation.enabled) {
        INPUT_VALIDATION(population, files, refpanel)

        QUALITY_CONTROL(
            INPUT_VALIDATION.out.validated_files,
            INPUT_VALIDATION.out.validation_report,
            site_files_ch.collect(),
            refpanel,
            params.build,
            population,
        )

        // check if QC chunks exist in case QC failed
        QUALITY_CONTROL.out.qc_metafiles.ifEmpty {
            error("QC step failed")
        }

        if (params.mode == 'imputation') {
            phased_ch = QUALITY_CONTROL.out.qc_metafiles

            if (phasing_engine != 'no_phasing') {
                PHASING(QUALITY_CONTROL.out.qc_metafiles, refpanel)
                phased_ch = PHASING.out.phased_ch
            }

            IMPUTATION(refpanel, params.imputation, phased_ch)

            if (params.merge_results == true) {
                ENCRYPTION(IMPUTATION.out.groupTuple(), encryption_password, refpanel)
            }
        }
    }

    ancestry_enabled = (params.ancestry != null) && params.ancestry.enabled

    if (ancestry_enabled) {
        ANCESTRY_ESTIMATION(files, params.ancestry, params.build)
    }

    if (params.pgs.enabled) {
        PGS_CALCULATION(
            IMPUTATION.out,
            ancestry_enabled ? ANCESTRY_ESTIMATION.out : channel.empty(),
        )
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

        println("::error:: Imputation job ${statusMessage}.")
        exit()
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
            println("::message:: Data have been exported successfully. We have sent a notification email to <b>${params.user.email}</b>")
        }
        else {
            println("::message:: Data have been exported successfully. We encrypted the results with the following password: <b>${encryption_password}</b>")
        }
        exit()
    }

    //PGS job
    if (params.send_mail && params.user.email != null) {
        sendMail(
            to: "${params.user.email}",
            subject: "[${params.service.name}] Job ${params.project} is complete",
            body: "Dear ${params.user.name},\n\nYour PGS job has finished successfully.\n\nYou can download the results from the following link: ${params.service.url}/index.html#!jobs/${params.project_id}",
        )
        println("::message:: Data have been exported successfully. We have sent a notification email to <b>${params.user.email}</b>")
    }
    else {
        println("::message:: Data have been exported successfully.")
    }
}
