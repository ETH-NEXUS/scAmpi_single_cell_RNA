
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


.. contents:: Table of Contents
   :depth: 2

Metacell Purity Analysis Report
===============================

.. image:: results/metacells/seacells_cells_per_metacell.png
.. image:: results/metacells/seacells_celltype_purity.png
{% set figure_counter = 1 %}
{% for sample in snakemake.config["samples"] %}

{{sample}}
-------------

.. list-table::
   :width: 100%
   :class: borderless

   * - .. image:: results/plotting/{{sample}}_seacells.first_celltype.png
          :width: 100%         
     - .. image:: results/plotting/{{sample}}.first_celltype.png
          :width: 100%

   * - .. image:: results/plotting/{{sample}}_seacells.phenograph.png
          :width: 100%         
     - .. image:: results/plotting/{{sample}}.phenograph.png
          :width: 100%

.. figure:: results/metacells/{{sample}}_seacells_celltype_hist.png
   :width: 90%
   :align: center    

**Figure {{ loop.index }}:** Metacell plots {{figure_counter}} for {{sample}}.   

.. include:: workflow/report/rules/seacells/{{sample}}_celltyping_summary.rst


{% set figure_counter = figure_counter + 1 %}
{% endfor %}

References
==========

.. [scAmpi2022] Bertolini A, et al. (2022) 
   scAmpi—A versatile pipeline for single-cell RNA-seq analysis from basics to clinics. PLoS Comput Biol 18(6): e1010097.
   (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010097).

