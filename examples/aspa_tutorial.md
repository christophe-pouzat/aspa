---
title: _Analyse des séquences de potentiels d'action_ tutorial
author: Christophe Pouzat
institute: Université Paris-Descartes and Centre national de la recherche scientitfique
documentclass: scrartcl
fontfamily: fourier
linkcolor: orange
---

# Data used

We are going to use spike trains obtained from the antennal lobe--first olfactory relay--of locust, _Schistocerca americana_. These spike trains can be found on the [zenodo-locust-datasets-analysis](https://christophe-pouzat.github.io/zenodo-locust-datasets-analysis/) GitHub repository. You can also find there a complete description of the sorting procedure used to go from the raw data, that are available on [zenodo](https://zenodo.org/record/21589), to the spike trains. We will mostly use the spike trains from experiment `locust20010214` that can be found at the following address: <https://github.com/christophe-pouzat/zenodo-locust-datasets-analysis/tree/master/Locust_Analysis_with_R/locust20010214/locust20010214_spike_trains>.

## Getting a spike train

We will start by downloading the spike train from unit 1 from `group` `Spontaneous_1`. This is done by typing in the shell:

~~~{#download-locust20010214_Spontaneous_1_tetB_u1 .bash}
wget https://raw.githubusercontent.com/christophe-pouzat/\
zenodo-locust-datasets-analysis/master/Locust_Analysis_with_R/\
locust20010214/locust20010214_spike_trains/\
locust20010214_Spontaneous_1_tetB_u1.txt
~~~
