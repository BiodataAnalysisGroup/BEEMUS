#!/usr/bin/env nextflow

// Script parameters
params.dataFolder = "data/test_new_data"
params.refFolder = params.dataFolder + "/ref"
params.vcfsFolder = params.dataFolder + "/vcf-files"
params.lineagesFolder = params.dataFolder + "/lineages"

refFolder = file(params.refFolder, type: 'dir')
refFile = file(params.refFolder + "/*.fasta")[0]
dataFolder = file(params.dataFolder, type: 'dir')
vcfsFolder = file(params.vcfsFolder, type: 'dir')
lineagesFolder = file(params.lineagesFolder, type: 'dir')
lineagesClasificationFile = file(params.dataFolder + "/SARS-CoV-2 lineage meta data.csv")
genesCoordinatesFile = file(params.refFolder + "/NC_045512.2_annot.gff3")

datasets = Channel
                .fromPath(params.vcfsFolder + '/*ann.vcf')
                .map { file -> tuple(file.baseName, file) }
datasets.into{
    dataset1
    dataset2
}


// process vcfs_preparation {
//     conda 'conda_envs/samtools_environment.yml'
//     publishDir 'data/test_new_data/indices', mode: 'copy'

//     input:
//     set datasetID, path(datasetFile) from dataset1

//     output:
//     path("${datasetID}.vcf.gz") into compressed_vcfs
//     path("${datasetID}.vcf.gz.csi") into compressed_vcfs_csi
    
    
//     script:
//     """
//     bgzip ${datasetFile}
//     tabix -f ${datasetID}.vcf.gz --csi > ${datasetID}.vcf.gz.csi
//     """
// }

// process AlternateReference {
//     conda 'conda_envs/gatk_environment.yml'
    

//     input:
//     set datasetID, path(datasetFile) from dataset2

//     output:
//     set datasetID, path("${datasetID}.fasta") into alternateFastaFiles

//     """
//     gatk IndexFeatureFile -I ${datasetFile} --verbosity ERROR
//     gatk FastaAlternateReferenceMaker -R $refFile -V ${datasetFile} -O ${datasetID}.fasta --verbosity ERROR
//     """
// }

// process renameHeadersFastaFiles {
//     publishDir 'data/test_new_data/AlternateReference', mode: 'copy'

//     input:
//     set datasetID, path(altfasta) from alternateFastaFiles

//     output:
//     path("${datasetID}.fasta") into alternateFastaFilesModified

//     """
//     tmp=${datasetID}
//     awk -v var="\${tmp%%_*}" -i inplace '/^>/{print ">" ++i "_" var; next}{print}' ${altfasta}
//     """
// }

// process combineFastaFiles {
//     publishDir 'data/test_new_data', mode: 'copy'

//     input:
//     path(altfasta_mod) from alternateFastaFilesModified.collect()
    
//     output:
//     path("combined.fasta") into combinedFastaFile

//     """
//     cat $altfasta_mod > combined.fasta
//     """
// }

process data_pre_processing {
    conda 'conda_envs/sars-cov-2_environment.yml'

    output:
    stdout into result

    script:
    """
    preprocessing.py --datafolder $dataFolder --vcfs $vcfsFolder --lineagesfiles $lineagesFolder --lineagesClasification "$lineagesClasificationFile" --genes_coordinates_path $genesCoordinatesFile
    """
}

// result.subscribe { println it }
