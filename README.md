# tasmitsThesis

<!-- ![logo](https://) -->

Quick overview
==============
This is a private repository containing Anastasios Mitsigkolas master thesis/internship

## Contents
- [Rationale](#rationale)
- [Quick start](#quick-start)
- [Contributing](#contributing)
- [License](#license)
- [Citations](#citations)
- [Credits](#credits)


Rationale
=========





Quick start
============

1. Conda or Miniconda package manager installed is required on your system.

2. Nextflow version 0.30.x or higher is required for the installation procedure, just download the distribution package by copying and pasting this command in your terminal:

    ```
    curl -fsSL https://get.nextflow.io | bash
    ```
    
    It creates the ``nextflow`` executable file in the current directory. You may want to move it to a folder accessible from your ``$PATH``.


3. Clone this GitHub repository:

    ```
    git clone https://github.com/
    ```

4. Create the folders structure inside the project's directory as specified below with the exact same names:

        .
        ├── ...
        ├── project's directory 
        │   ├── data                             # Data folder
        │   │   ├── vcf-files                    # The folder where all vcf files live in
        │   │   ├── ref                          # The folder where the reference genome and its indices live in.
        │   │   │   ├── NC_045512.2.fasta
        │   │   │   ├── NC_045512.2.fasta.fai
        │   │   │   ├── NC_045512.2.dict
        │   │   │   ├── NC_045512.2_annot.gff3
        │   │   │   └── ...
        │   │   ├── biosample_result.csv         # The metadata file
        │   │   └── lineage-classification.csv   # The samples to lineages classification file
        │   └── ...
        └── ...
5. Create and activate conda environment:

    ```
    conda env create -f conda_envs/sars-cov-2_environment.yml
    ```

6. Run the pipeline:

    ```
    nextflow pipeline.nf
    ```

7. Inside the data folder you will find a summary file ``dataset_summary.csv``. Select the most frequently occurring lineages along your dataset and write them down for the next step.

8. Activate the conda environment:

    ```
    conda activate sars_work
    ```

9. Run the graph analysis script (from the project's folder), replacing the string after the ``--lineages`` flag, with the lineages of interest in sequential order, separated by space, as shown below: 

    ``--lineages B.1.1.7 B.1.617.2 AY.43 BA.1 AY.12 AY.4 AY.122 AY.9 BA.1.1``
    ```
    python bin/graph_analysis.py \
        --datafolder data/ \
        --datasetname dataset.csv \
        --graphsfolder data/graphs/ \
        --lineages <lineages in sequential order separated by space> \
        --export_arrow_edges \
        --export_nodes_long_format
    ```

10. If everything went well, a browser window should automatically pop up, showing the produced graph, like this depicted below, in a non-arranged form, that you should rearrange by drag and drop the nodes based on your preferences. After that we strongly recommend writing down the unique ids of the nodes of interest, separated by space (shown underneath the node's lineage(s)), in groups of interest:

    ``E.g. <11 7 9 8> <51 90 72>``

    ![Alt text](imgs/graph.png?raw=true "Graph")

11. Run the msa analysis script (from the project's folder), replacing the string after the ``--nodes_of_interest`` ``--lineages`` flag, with the lineages and the nodes of interest in sequential order, separated by space, as shown below:

    ``--lineages B.1.1.7 B.1.617.2 AY.43 BA.1 AY.12 AY.4 AY.122 AY.9 BA.1.1``
    ``--nodes_of_interest 39 34 86``    
    ```
    python bin/msa.py \
        --datafolder data \
        --datasetname dataset.csv \
        --msafilename MSA.fasta \
        --graphsfolder data/graphs \
        --nodesfilename nodes \
        --treefilepath ../tree-of-life/files/life.json \
        --lineages <lineages in sequential order separated by space> \
        --nodes_of_interest <nodes of interest in sequential order separated by space>
    ```

    - You can run the ``msa.py`` script as many times as it is necessary, while using different nodes of interest.
    - Caution! Remember to change the ``--treefilepath`` file name along different running, since it will replace the previous file each time it executes.

12. Go to the ``tree-of-life`` folder (from the project's folder) and start a local http server using python:

    ```
    cd ../tree-of-life
    python -m http.server
    ```
13. Open a browser of your choice and type ``http://127.0.0.1:8000/``. If everything is all right, you will be able to see a tree graph such as the following:

    ![Alt text](imgs/tree-example.png?raw=true "Tree of Life")

Contributing
============

Project contribution are more than welcome. 


License
=======

This project is released under the MIT license. See the [LICENSE](LICENSE.md) file for details.

Citations
=========

If you use this project in your research, please cite:


Credits
=======