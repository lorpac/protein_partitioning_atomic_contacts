# Partitioning of atomic contacts in a protein variants

Data and code for the comparison the local partitioning of atomic contacts in the Amino Acid Networks of protein variants. The presented method and results are described in:
- Lorenza Pacini, Claire Lesieur, A computational methodology to diagnose sequence-variant dynamic perturbations by comparing atomic protein structures, Bioinformatics, Volume 38, Issue 3, 1 February 2022, Pages 703–709, https://doi.org/10.1093/bioinformatics/btab736

## Requirements

Python 3.x is required and [jupyter](https://jupyter.org/) is needed to run the `case_studies` notebook. The calculations are built upon Rodrigo Dorantes-Gilardi's module [biographs](https://github.com/rodogi/biographs). The modules required for the calculation and the plotting of the results are listed in `requirements.txt` and can be installed using

```
pip install -r requirements.txt
```
You need to have LaTex installed on your system in order to produce the plots with Matplotlib. On Ubuntu, you can install LaTex and the necessary extensions by running

```
sudo apt-get install dvipng texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```
and if you still get an error, try installing `cm-super`:
```
sudo apt-get install cm-super
```

For other operating systems, or if you encounter problems, please follow the instructions in Matplotlib's [tutorial](https://matplotlib.org/3.1.0/tutorials/text/usetex.html).


## Author

Lorenza Pacini - [lorpac](https://github.com/lorpac)

## How to cite

If you use this code in your work, please cite:
- Lorenza Pacini, Claire Lesieur, A computational methodology to diagnose sequence-variant dynamic perturbations by comparing atomic protein structures, Bioinformatics, Volume 38, Issue 3, 1 February 2022, Pages 703–709, https://doi.org/10.1093/bioinformatics/btab736
  
## Licence

The Building network source code is available under the [CeCILL](http://cecill.info/) licence. Please see `LICENCE.txt` for details.

