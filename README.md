# SENMAP
SENMAP: A Convolutional Neural Network Architecture for Curation of LTR-RT Libraries from Plant Genomes
Transposable elements (TEs) are specific structures of the genome of species, which can move from one location to another. For that reason, they can cause mutations or changes that can be negative, such as the appearance of diseases, or beneficial, such as participating in fundamental roles in the evolution of genomes and  genetic diversity. Long Terminal Repeat retrotransposons (LTR-RT) are the most abundant in plant species, hence the importance of studying these structures in particular. Over the time, these elements can suffer changes called nested insertions, which can inactivate or modify the functioning of the element, for that they are no longer consider as intact element and cannot be used for identification and classification studies. In this work we present SENMAP, a convolutional neural network architecture to obtain intact LTR-RT sequences in plant genomes, which is composed by four convolutional layers, LeakyReLU as activation function and  BinaryFocalLoss as loss function. Achieving an F1-score percentage of 91.37% with test data, identifying low quality sequences rapidly and efficiently, contributing to curate libraries of LTR retrotransposons of plants genomes published in large-scale sequencing projects due to the post-genomic era. 

## Usage
First create the output folder, and then execute SENMAP:

```
mkdir output_folder
```

```
python3 SENMAP.py -l LTR_sequences_lib.fa -o output_folder
```

## Citation
if you use this software, the neural network architecture, or any other part of the code, please cite us as following:

* Orozco-Arias, S., Candamil-Cortés, M. S., Valencia-Castrillón, E., Jaimes, P. A., Orozco, N. T., Arias-Mendoza, M., ... & Isaza, G. (2021, October). SENMAP: A Convolutional Neural Network Architecture for Curation of LTR-RT Libraries from Plant Genomes. In 2021 IEEE 2nd International Congress of Biomedical Engineering and Bioengineering (CI-IB&BI) (pp. 1-4). IEEE.

