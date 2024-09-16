
General overview
----------------

This scAmpi workflow is organized into two main parts: 
the ``scAmpi_basic`` part and the ``scAmpi_clinical`` part, 
which can be run independently. 
scAmpi_basic includes general scRNA processing steps, 
such as mapping, QC, normalisation, unsupervised clustering, 
cell type classification, and DE analysis. 
For more details see the  publication
[scAmpi2022]_.

scAmpi_clinial includes the search for disease relevant drug targets for differentially expressed genes. 
Note that the clinical part is only applied if at least one cluster identified in your sample is indicated 
as a diseased ("malignant") cell type.


.. [scAmpi2022] Bertolini A, et al. (2022) 
   scAmpi—A versatile pipeline for single-cell RNA-seq analysis from basics to clinics. PLoS Comput Biol 18(6): e1010097.
   (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010097).

Metacell Purity Analysis Report
===============================
{% for sample in snakemake.config["samples"] %}

{{sample}}
-------------

.. include:: workflow/report/rules/seacells/{{sample}}_celltyping_summary.rst
   
.. image:: results/metacells/{{sample}}_seacells_celltype_dens.png
    :width: 600
{% endfor %}