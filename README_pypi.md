Scoary - Microbial pan-GWAS
###########################

`Scoary <https://github.com/AdmiralenOla/Scoary>`_ is designed to take the gene_presence_absence.csv file from `Roary <https://sanger-pathogens.github.io/Roary/>`_ as well as a traits file created by the user and calculate the assocations between all genes in the accessory genome and the traits. It reports a list of genes sorted by strength of association per trait.

Dependencies
------------
- Python (Tested with versions 2.7, 3.4, 3.5, 3.6 and 3.6-dev)
- `SciPy <http://www.scipy.org/install.html>`_ (Tested with versions 0.16, 0.17, 0.18)

If you supply custom trees (Optional)

- ete3
- six

Note that ete3 and six are not automatically installed. You can do `pip install ete3 six` to get them

Using the GUI (Optional)

- Tkinter/ttk

Tkinter/ttk is already part of most python distributions. If you lack it consider getting Homebrew/Linuxbrew and running `brew install python --with-tcl-tk`

Installation
------------
The easiest way to install Scoary is through the pip package manager::

    pip install scoary

Usage
-----
    scoary -g <gene_presence_absence.csv> -t <traits.csv>

Documentation
-------------
The most updated documentation for scoary is found at `the project site <https://github.com/AdmiralenOla/Scoary>`_

Citation
--------
If you use Scoary, please cite our `paper <https://dx.doi.org/10.1186/s13059-016-1108-8>`_

License
-------
Scoary is freely available under a GPLv3 license.

Contact
-------
Ola Brynildsrud (ola.brynildsrud@fhi.no)
