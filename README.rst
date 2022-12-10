|Stars| |PyPI| |Docs| |Total downloads| |Monthly downloads|

.. |Stars| image:: https://img.shields.io/github/stars/genecell/COSG?logo=GitHub&color=yellow
   :target: https://github.com/genecell/COSG/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/cosg?logo=PyPI
   :target: https://pypi.org/project/cosg
.. |Docs| image:: https://readthedocs.org/projects/cosg/badge/?version=latest
   :target: https://cosg.readthedocs.io
.. |Total downloads| image:: https://static.pepy.tech/personalized-badge/cosg?period=total&units=international_system&left_color=black&right_color=orange&left_text=downloads
   :target: https://pepy.tech/project/cosg
.. |Monthly downloads| image:: https://static.pepy.tech/personalized-badge/cosg?period=month&units=international_system&left_color=black&right_color=orange&left_text=downloads/month
 :target: https://pepy.tech/project/cosg

Accurate and fast cell marker gene identification with COSG
=======================================================================================================

Overview
---------
COSG is a cosine similarity-based method for more accurate and scalable marker gene identification.

- COSG is a general method for cell marker gene identification across different data modalities, e.g., scRNA-seq, scATAC-seq and spatially resolved transcriptome data.
- Marker genes or genomic regions identified by COSG are more indicative and with greater cell-type specificity.
- COSG is ultrafast for large-scale datasets, and is capable of identifying marker genes for one million cells in less than two minutes.

The method and benchmarking results are described in `Dai et al., (2021)`_. 

Documentation
--------------
The documentation for COSG is available `here <https://cosg.readthedocs.io/en/latest/>`_.

Tutorial
---------

The `COSG tutorial <https://nbviewer.jupyter.org/github/genecell/COSG/blob/main/tutorials/COSG-tutorial.ipynb>`_ provides a quick-start guide for using COSG and demonstrates the superior performance of COSG as compared with other methods, and the `Jupyter notebook <https://github.com/genecell/COSG/blob/main/tutorials/COSG-tutorial.ipynb>`_ is also available.

Question
---------
For questions about the code and tutorial, please contact Min Dai, daimin@zju.edu.cn.

Citation
---------
If COSG is useful for your research, please consider citing `Dai et al., (2021)`_.

.. _Dai et al., (2021): https://www.biorxiv.org/content/10.1101/2021.06.15.448484v1


