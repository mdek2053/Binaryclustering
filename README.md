# Binary Clustering

This repository contains code related to a paper for my bachelor thesis([link](https://repository.tudelft.nl/islandora/object/uuid%3A9c11c1a3-fe0d-4be9-89dc-7519dba3e7de?collection=education)).

The main idea is to apply a clustering algorithm to binary vectors instead of integer vectors, in order to achieve a significant speedup. This is because two binary vectors can be compared very quickly using native CPU instructions like XOR. 

# Contributors
The algorithm was mostly developed in a team of 4:

- Milan de Koning
- Jurriën Theunisz
- Bartek Golik
- Pavel Verigo

The integration with R and experiments were created and conducted by me.