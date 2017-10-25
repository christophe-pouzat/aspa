<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#the-idea-motivating-this-development">1. The idea motivating this development</a></li>
<li><a href="#required-software">2. Required software</a>
<ul>
<li><a href="#getting-and-compiling-the-code">2.1. Getting and compiling the code</a></li>
</ul>
</li>
<li><a href="#data-used">3. Data used</a>
<ul>
<li><a href="#getting-a-spike-train">3.1. Getting a spike train</a></li>
</ul>
</li>
<li><a href="#preliminary-analysis">4. Preliminary analysis</a>
<ul>
<li><a href="#reading-the-data">4.1. Reading the data</a></li>
<li><a href="#basic-statistics">4.2. Basic statistics</a></li>
<li><a href="#basic-plots">4.3. Basic plots</a>
<ul>
<li><a href="#org7e606ea">4.3.1. Raster plot</a></li>
<li><a href="#org0174e70">4.3.2. A fancy trick</a></li>
<li><a href="#org3574d03">4.3.3. Counting process plots</a></li>
<li><a href="#orga5646da">4.3.4. Lagged ranked ISI plot</a></li>
</ul>
</li>
<li><a href="#org890c8ed">4.4. Histograms</a>
<ul>
<li><a href="#org2506a23">4.4.1. Getting the cross validation score</a></li>
<li><a href="#org26881d6">4.4.2. Building the histogram</a></li>
</ul>
</li>
<li><a href="#orgd170863">4.5. ISI histograms of the seven "good" neuron of dataset <code>locust20010214</code></a></li>
</ul>
</li>
</ul>
</div>
</div>


<a id="the-idea-motivating-this-development"></a>

# The idea motivating this development

The idea here is to implement the Unix/Linux "philosophy"&#x2013;as exposed
for instance in the article of Arnold Robbins
[*What's GNU*](http://www.linuxjournal.com/article/2762)&#x2013;to the
analysis of neuronal spike trains. Since spike trains make not too
voluminous data, they can be stored as text files (ASCII) and most
operations on them can be designed as "filters", that is programs
(usually but not always written in `C`) that read their input in text
format from the "standard input" (`stdin`) and send their result in text
format to the "standard output" (`stdout`). For graphical displays, we
are going to use [gnuplot](http://gnuplot.info/).


<a id="required-software"></a>

# Required software

The code will be written mostly in `C`. If you want a clear and quick
introduction to this language, check Ben Klemens:
[Modeling With Data](http://modelingwithdata.org/about_the_book.html).
To compile the code you will need a `C` compiler like
[gcc](https://gcc.gnu.org/). If you are using `Linux` or `MacOS` it's
in a package from your favorite distribution, if you are using `Windows`
you will have to install [Cygwin](https://cygwin.com/index.html). The
heavy computational work is going to be performed mainly by the
[gsl](http://www.gnu.org/software/gsl/) (the *GNU Scientific Library*)
that is easily installed through your package manager (from now on, for
windows users, the "package manager" refers to the one of `Cygwin`). The
graphs will be generated with [gnuplot](http://www.gnuplot.info/).
Windows user who want to use the interactive plotting capabilities of
the library (recommended) will also have to install
[`Cygwin/X`](http://x.cygwin.com/).

For now the compilation requires either `make` or [Scons](http://scons.org/).


<a id="getting-and-compiling-the-code"></a>

## Getting and compiling the code

The code is hosted on
[`GitHub`](https://github.com/christophe-pouzat/aspa). The easiest is
to clone or download the repository (there is a button for that on the
GitHub page). Once you have the repository on your hard drive, go to the
code sub-directory and, if using `make`, type:

    cd ../code && make all

    cc `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11    -c -o aspa_single.o aspa_single.c
    ar cr libaspa.a aspa_single.o
    cc `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11    -c -o aspa_read_spike_train.o aspa_read_spike_train.c
    cc aspa_read_spike_train.o libaspa.a `pkg-config --libs gsl `  -o aspa_read_spike_train
    cc `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11    -c -o aspa_mst_fns.o aspa_mst_fns.c
    cc aspa_mst_fns.o libaspa.a `pkg-config --libs gsl `  -o aspa_mst_fns
    cc `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11    -c -o aspa_mst_aggregate.o aspa_mst_aggregate.c
    cc aspa_mst_aggregate.o libaspa.a `pkg-config --libs gsl `  -o aspa_mst_aggregate
    cc `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11    -c -o aspa_mst_plot.o aspa_mst_plot.c
    cc aspa_mst_plot.o libaspa.a `pkg-config --libs gsl `  -o aspa_mst_plot

or with `SCons`:

    scons -Q

This will compile the library `libaspa.a` as well as a bunch of user
programs all starting with `aspa_`, like `aspa_read_spike_train`. As
mentioned previously, you need the `gsl` to be installed in order to
compile the code.

Once the compilation is done you should move the user programs to one of
the directories listed on your `PATH`, that is on one of the directories
appearing when you type:

    echo $PATH

    /opt/jython/bin/:/usr/local/sbin:/usr/local/bin:/usr/bin:/usr/lib/jvm/default/bin:/usr/bin/site_perl:/usr/bin/vendor_perl:/usr/bin/core_perl

After that, you're in business.


<a id="data-used"></a>

# Data used

We are going to use spike trains obtained from the antennal lobe&#x2013;the first
olfactory relay&#x2013;of locusts, *Schistocerca americana*. These spike trains
can be found on the
[zenodo-locust-datasets-analysis](https://christophe-pouzat.github.io/zenodo-locust-datasets-analysis/)
GitHub repository. You can also find there a complete description of the
sorting procedure used to go from the raw data, that are available on
[zenodo](https://zenodo.org/record/21589), to the spike trains. We
will mostly use the spike trains from experiment `locust20010214` that
can be found at the following address:
<https://github.com/christophe-pouzat/zenodo-locust-datasets-analysis/tree/master/Locust_Analysis_with_R/locust20010214/locust20010214_spike_trains>.


<a id="getting-a-spike-train"></a>

## Getting a spike train

We will start by downloading the spike train from unit 1 from `group`
`Spontaneous_1`. This is done by typing in the shell (I'm using the
"line continuation character, " to fit my lines on a single page of the
`PDF` version of this document, when typing directly to the shell you
don't need these line breaks):

    wget https://raw.githubusercontent.com/christophe-pouzat/\
    zenodo-locust-datasets-analysis/master/Locust_Analysis_with_R/\
    locust20010214/locust20010214_spike_trains/\
    locust20010214_Spontaneous_1_tetB_u1.txt

    --2017-10-23 21:14:47--  https://raw.githubusercontent.com/christophe-pouzat/zenodo-locust-datasets-analysis/master/Locust_Analysis_with_R/locust20010214/locust20010214_spike_trains/locust20010214_Spontaneous_1_tetB_u1.txt
    Certificat de l'autorité de certification « /etc/ssl/certs/ca-certificates.crt » chargé
    Résolution de raw.githubusercontent.com… 151.101.92.133
    Connexion à raw.githubusercontent.com|151.101.92.133|:443… connecté.
    requête HTTP transmise, en attente de la réponse… 200 OK
    Taille : 27743 (27K) [text/plain]
    Sauvegarde en : « locust20010214_Spontaneous_1_tetB_u1.txt »
    
         0K .......... .......... .......                         100%  957K=0,03s
    
    2017-10-23 21:14:47 (957 KB/s) — « locust20010214_Spontaneous_1_tetB_u1.txt » sauvegardé [27743/27743]

This "spike train" contains in fact the result of 30 consecutive
continuous acquisitions, each 29 s long with a 1 s gap in between, as is
made clear in the
[detailed
sorting description](https://christophe-pouzat.github.io/zenodo-locust-datasets-analysis/Locust_Analysis_with_R/locust20010214/Sorting_20010214_tetB.html) of this data set.


<a id="preliminary-analysis"></a>

# Preliminary analysis


<a id="reading-the-data"></a>

## Reading the data

In is not expected that the data (spike trains) one wants to work with
will be obtained in any standard format. That means that a usually
slightly "painful" work will be required (but that's always the case
when dealing with actual data) to read the data and reformat them in the
text (or binary) format used by `aspa`. Looking at the source code of
`aspa_read_spike_train` is the way to proceed (more specifically, look
at the code of `aspa_raw_fscanf` that is called by
`aspa_read_spike_train` and that is found in `aspa_single.c`).

The data we just downloaded are collections of spike times in "sample
times"&#x2013;the time unit is therefore not the second but 1/15000
second&#x2013;with one spike time per line. This can be seen by calling first
the `head` function (showing by default the first ten lines of the
file):

    head locust20010214_Spontaneous_1_tetB_u1.txt

    4364.629
    49876.8
    50529.95
    50988.26
    51371.66
    51769.29
    52703.77
    54772.34
    56472.7
    71766.51

Calling `tail` shows the last lines of the file (by default the last ten
lines):

    tail locust20010214_Spontaneous_1_tetB_u1.txt

    13442792
    13455679
    13458610
    13460049
    13460517
    13461154
    13464139
    13470059
    13471539
    13472243

Function `aspa_read_spike_train` will read these times from the `stdin`
and output them in a "nice" format (still a text file by default) to the
`stdout`. You can get a description to arguments accepted by the
function by calling it with the `--help` argument:

    ./aspa_read_spike_train --help

    Usage: 
      --in_bin: specify binary data input
      --out_bin: specify binary data output
      --sample2second <positive real>: the factor by which times
      in input data are divided in order get spike times in seconds
      used only when reading 'raw' data (default 15000)
      --inter_trial_interval <positive real>: the inter trial
      interval (in s) used only when reading 'raw' data
      --trial_duration <positive real>: the recorded duration
      (in s) of each trial used only when reading 'raw' data
      --stim_onset <real>: the stimulus onset time
      (in s) if that makes sense, used only when reading 'raw' data
      --stim_offset <real>: the stimulus offset time
      (in s) if that makes sense, used only when reading 'raw' data

For demonstration we can call it on the data file we just downloaded
(`locust20010214_Spontaneous_1_tetB_u1.txt`), writing the result into a
new text file `locust20010214_Spontaneous_1_tetB_u1.aspa` for further
inspection:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 < \
    locust20010214_Spontaneous_1_tetB_u1.txt > \
    locust20010214_Spontaneous_1_tetB_u1.aspa

We can then look at the first 25 lines of our new file with:

    head -n 25 locust20010214_Spontaneous_1_tetB_u1.aspa 

    # Number of trials: 28
    # Number of aggregated trials: 1
    # Stimulus onset: 0 (s)
    # Stimulus offset: 0 (s)
    # Single trial duration: 29 (s)
    
    
    # Start of trial: 0
    # Trial start time: 0 (s)
    # Number of spikes: 94
    0.290975
    3.32512
    3.36866
    3.39922
    3.42478
    3.45129
    3.51358
    3.65149
    3.76485
    4.78443
    5.06381
    5.11507
    5.24077
    5.28448
    5.31933

We see that the "non-data" element are on lines starting with a "#"
character. The "head" of the file specifies how many trial are in the
file and gives some other information. The data from trial 0 (we start
counting at 0) com next after two blank lines. To see the whole file
interactively you can type:

    less locust20010214_Spontaneous_1_tetB_u1.aspa 


<a id="basic-statistics"></a>

## Basic statistics

Program `aspa_mst_fns` (`mst` stands for "multiple spike trains" and
`fns` for "[Five-number summary](https://en.wikipedia.org/wiki/Five-number_summary)") return elementary statics related to a spike train data set.
A description of its use is obtained by calling the program with the `--help` argument:

    ./aspa_mst_fns --help

    Usage: 
      --in_bin: specify binary data input
    
    Returns five number summary and additional stats.

We can call this function directly on the output of `aspa_read_spike_train` using a [pipe](http://www.linfo.org/pipe.html) with:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 < \
    locust20010214_Spontaneous_1_tetB_u1.txt | \
    ./aspa_mst_fns

    Data from 28 trials.
    The mean rate is: 4.10222 Hz.
    The inter spike interval statistics are:
      The sample contains 3303 elements.
      The mean and SD are   : 0.2333 and 0.4660.
      The median and MAD are: 0.0546 and 0.0359.
    The five number summary:
      Min.   1st qrt Median 3rd qrt Max. 
      0.0157 0.0369  0.0546 0.1491  4.5264
    A 95% confidence interval for the lag 1 Spearman rank correlation is: [0.400336,0.443483].


<a id="basic-plots"></a>

## Basic plots

There are several plots one might want to create at an early stage of spike train data analysis. Since most of these plots are more "attractive" when built from data with a response to a stimulus, we will start by getting one such case (from the same experiment and same neuron):

    wget https://raw.githubusercontent.com/christophe-pouzat/\
    zenodo-locust-datasets-analysis/master/Locust_Analysis_with_R/\
    locust20010214/locust20010214_spike_trains/\
    locust20010214_C3H_1_tetB_u1.txt

    --2017-10-23 20:53:00--  https://raw.githubusercontent.com/christophe-pouzat/zenodo-locust-datasets-analysis/master/Locust_Analysis_with_R/locust20010214/locust20010214_spike_trains/locust20010214_C3H_1_tetB_u1.txt
    Certificat de l'autorité de certification « /etc/ssl/certs/ca-certificates.crt » chargé
    Résolution de raw.githubusercontent.com… 151.101.92.133
    Connexion à raw.githubusercontent.com|151.101.92.133|:443… connecté.
    requête HTTP transmise, en attente de la réponse… 200 OK
    Taille : 29318 (29K) [text/plain]
    Sauvegarde en : « locust20010214_C3H_1_tetB_u1.txt.1 »
    
         0K .......... .......... ........                        100% 1,75M=0,02s
    
    2017-10-23 20:53:01 (1,75 MB/s) — « locust20010214_C3H_1_tetB_u1.txt.1 » sauvegardé [29318/29318]

This file contains the responses to 25 stimulations with `cis-3-hexen-1-ol`. The classical way of displaying such data is the `raster plot`. This plot as well as several over ones we will shortly see is generated by calling `aspa_mst_plot`. As usual, calling the function with argument `--help` gives us a basic explanation on how to use it:

    ./aspa_mst_plot --help

    Usage: 
      --in_bin: specify binary data input
      --text: specify text output
      --what <string>: one of 'raster', 'cp_rt', 'cp_wt',
      'cp_norm', 'lrank', the type of plot (see bellow)
      --lag <positive integer>: the lag used in lagged
        ranked plots (default at 1).
    
    An interactive plot is generated.
    If what is set to 'raster' a raster plot is generated.
    If what is set to 'cp_rt' the observed counting process
    in 'real' time is generated, that is trial appear one after
    the other.
    If what is set to 'cp_wt' the observed counting processes
    corresponding to each trial are shown on the 'within trial time'.
    If what is set to 'cp_norm' the normalized aggregated counting
    process is displayed (normalization means here that the step size
    due to each spike in each trial is 1/number of trials; in a sense
    the 'mean' counting process is displayed).
    If what is set to 'lrank', isi are ranked from the smallest to
    the largest and the rank of isi i+lag is plotted against the
    lag of isi i.


<a id="org7e606ea"></a>

### Raster plot

Here, to get the classical raster we do:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 \
    --stim_onset=10 --stim_offset=11 < \
    locust20010214_C3H_1_tetB_u1.txt | \
    ./aspa_mst_plot --what=raster

This will make a new window appear with a plot similar to the one we will now construct after calling the function with an additional argument (you can type `q` to kill the plot window):

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 \
    --stim_onset=10 --stim_offset=11 < \
    locust20010214_C3H_1_tetB_u1.txt | \
    ./aspa_mst_plot --what=raster --text > \
    locust20010214_C3H_1_tetB_u1.raster

Here instead of the "new window output" we generated at text output (that's what the `--text` argument means) sent to the `stdout` and redirected this `stdout` to a file called `locust20010214_C3H_1_tetB_u1.raster`. We can now build "by hand" with `gnuplot` the same figure as the one we directly got (we have now more control on the output):

    set grid
    unset key
    set xlabel 'Time (s)'
    set ylabel 'Trial'
    plot [0:30] [0:26] 'locust20010214_C3H_1_tetB_u1.raster' \
         index 0 using 1:2 with filledcurve closed lc 'grey',\
         '' index 1:25 using 1:2 with dots lc 'black'

![img](fig/locust20010214_C3H_1_tetB_u1_raster.png)


<a id="org0174e70"></a>

### A fancy trick

We can also make the raster plot and get the basic stats printed at once with the [tee](https://www.gnu.org/software/coreutils/manual/html_node/tee-invocation.html) command as follows:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 < \
    locust20010214_C3H_1_tetB_u1.txt | tee >(aspa_mst_plot --what=raster) | \
    ./aspa_mst_fns


<a id="org3574d03"></a>

### Counting process plots

There are several ways to create a "counting process" plot. The first one, used mainly for checking data stationarity is building the "true" observed counting process plot, that is at each spike time the step function jumps by one unit and successive trials are shown one after the other as they *actually* occurred. This is what is specified with argument `cp_rt` to option `what`:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 \
    --stim_onset=10 --stim_offset=11 < \
    locust20010214_C3H_1_tetB_u1.txt | \
    ./aspa_mst_plot --what=cp_rt

Giving a plot looking like:

![img](fig/locust20010214_C3H_1_tetB_u1_cp_rt.png)

We might also want to look at the individual observed counting processes after realigning them on the stimulus onset. This is obtained with argument `cp_wt` to option `what`:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 < \
    locust20010214_C3H_1_tetB_u1.txt | \
    ./aspa_mst_plot --what=cp_wt

resulting in a plot looking like:

![img](fig/locust20010214_C3H_1_tetB_u1_cp_wt.png)

We can also decide that to see if there is a response or not, we can construct the average step function. That is, we replace the step size in the previous plot by 1/N (N is the number of trials) and we sum all these resulting step functions. This is done with argument `cp_norm` to option `what`:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 \
    --stim_onset=10 --stim_offset=11 < \
    locust20010214_C3H_1_tetB_u1.txt | \
    ./aspa_mst_plot --what=cp_norm

resulting in a plot looking like:

![img](fig/locust20010214_C3H_1_tetB_u1_cp_norm.png)

One might want to prepare a more sophisticated plot for, say, a publication showing both the individual observed counting processes and their average. This could be done as follows, by creating first some files with the "properly" formated data:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 \
    --stim_onset=10 --stim_offset=11 < \
    locust20010214_C3H_1_tetB_u1.txt > \
    locust20010214_C3H_1_tetB_u1.aspa
    
    ./aspa_mst_plot --what=cp_norm --text < \
    locust20010214_C3H_1_tetB_u1.aspa > \
    locust20010214_C3H_1_tetB_u1.cp_norm
    
    ./aspa_mst_plot --what=cp_wt --text < \
    locust20010214_C3H_1_tetB_u1.aspa > \
    locust20010214_C3H_1_tetB_u1.cp_wt

Then within `gnuplot`:

    set grid
    unset key
    set xlabel 'Time (s)'
    set title 'Observed counting processes'
    set ylabel 'Events count'
    plot [0:29] [0:245] 'locust20010214_C3H_1_tetB_u1.cp_wt' \
         index 0 using 1:2 with filledcurve closed lc 'grey',\
         '' index 1:25 using 1:2 with steps lc 'black',\
         'locust20010214_C3H_1_tetB_u1.cp_norm' index 1 \
         using 1:2 with steps lc 'red' linewidth 2

![img](fig/locust20010214_C3H_1_tetB_u1_cp_norm_wt.png)


<a id="orga5646da"></a>

### Lagged ranked ISI plot

Coming back to the spontaneous data, a good graphical way to look for correlations between successive inter spike intervals is to rank them (from the smallest to the largest) before plotting the rank of interval i+lag against the rank of interval i. We can do that with the spontaneous data as follows:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 --stim_onset=10 --stim_offset=11 < \ locust20010214_Spontaneous_1_tetB_u1.txt | ./aspa_mst_plot --what=lrank

We then obtain a graph looking like:

![img](fig/locust20010214_Spontaneous_1_tetB_u1_lrank.png)


<a id="org890c8ed"></a>

## Histograms

A cross-validation based binwidth selection is implemented in `aspa`. It follows Mats Rudemo "Empirical Choice of Histograms and Kernel Density Estimators" (*Scandinavian Journal of Statistics* **9**:65-78, 1982) approach. In short given a sample of inter spike intervals (ISI), a cross validation score approximating the integrated mean squared error (between the estimated density and the true one) is computed for a sequence of binwidths (or equivalently number of bins). This curve score vs number of bins can be plotted and the the "best" choice is the number of bins minimizing it. An example follows in the next sub-section.


<a id="org2506a23"></a>

### Getting the cross validation score

We first get the ISI sample using function `aspa_mst_isi` whose description is:

    ./aspa_mst_isi --help

    usage: ./aspa_mst_isi [-b --bin] [-h --help]
    
      -b --bin: the data read from the 'stdin' are in binary format.
      -h --help: prints this message.
     The program reads data from the 'stdin' (default in text format)
     most likely resulting from a call to 'aspa_read_spike_train',
     gets the inter spike intervals and writes the number of ISIs
     followed by the individual ISIs (one per line) to the stdout

We use it as follows, storing the obtained ISI sample in a file named `u1isi1.txt`:

    ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 < \
    			locust20010214_Spontaneous_1_tetB_u1.txt | \
        ./aspa_mst_isi > u1isi1.txt

We can check the head of the generated file with `head`:

    head u1isi1.txt

    3303
    3.03414
    0.04354
    0.03056
    0.0255599
    0.02651
    0.0622902
    0.13791
    0.11336
    1.01958

We use next function `aspa_hist_bw` whose description is:

    ./aspa_hist_bw --help

    usage: ./aspa_hist_bw [-l --log] [-f --from=integer]
              [-t --to=integer] [-b --best] [-h --help]
    
      -l --log: should the log of the observations be used?
      -f --from <positive integer>: the smallest number of bins to explore
         (default set to 2).
      -t --to <positive integer>: the largest number of bins to explore
         (default set to the number of observations divided by 5).
      -b --best: should only the best number of bins be printed to the 'stdout'?
      -h --help: prints this message.
     The program reads data from the 'stdin' (in text format),
     the first line should contain the number of observations (integer)
     the following lines should contain the observations, one per line
     in decimal notation. If a log transformed of the data is requested
     it is applied first. Then a sequence of histograms with number of
     bins variying between from and to (inclusive) is contructed. For each
     histogram, the empirical bin probability (bin count divided by sample
     size), p_hat_i for each bin i is obtained and the score is computed as
     follows: (2-(n+1) * sum(p_hat_i^2))*m/(n-1) where m is the number of
     bins and where the summation is done over the bins. The best number
     of bins is the one giving the smallest score.
     The program prints to the 'stdout' the number of bins and the
     corresponding cross-validation scores on two columns (default)
      or only the best number of bins if 'best' is used as an argument.

On the ISI sample we just got that gives (using the log of the observations to limit the consequences of the skewness of the distribution):

    ./aspa_hist_bw --log < u1isi1.txt > u1bw1.txt

    Sample size: 3303
    Exploring now number of bins between 2 and 330.
    Using a log transformation of the data.
    Histograms will be built between -4.15413 and 1.50993.
    The best number of bins is: 46 giving a score of -2.07509

Using `gnuplot` we quickly visualize the curve:

    set grid
    unset key
    set xlabel 'Number of bins'
    set ylabel 'Cross validation score'
    plot 'u1bw1.txt' using 1:2 with lines lc 'black'\
         linewidth 2

![img](fig/locust20010214_Spontaneous_1_tetB_u1_bw.png)


<a id="org26881d6"></a>

### Building the histogram

Once the "best" number of bins is known (or at least one has decided on the number of bins to use), we call `aspa_hist` to construct the histogram. Its description is:

    ./aspa_hist --help

    usage: ./aspa_hist -n --n_bins=integer [-l --log] 
              [-p --prob] [-h --help]
    
      -n --n_bins <positive integer>: the number of bins to use.
      -l --log: should the log of the observations be used?
      -p --prob: should the result be normalized so that the
         the histogram integral is one?
      -h --help: prints this message.
     The program reads data from the 'stdin' (in text format),
     the first line should contain the number of observations (integer)
     the following lines should contain the observations, one per line
     in decimal notation. If a log transformed of the data is requested
     it is applied first. Then a histogram with 'n_bins' bins is contructed.
     If a log transformation was used, the bins boundaries are transformed
     back to the original scale. If argument 'prob' is specified, the histogram
     is normalized such that its integral is one (it becomes of proper PDF
     estimator), otherwise the bin counts are kept.
     The program prints to the 'stdout' on 3 columns: the left bin boundary;
     the right bin boundary; the bin count or bin frequency (depending on the
     specification of argument 'prob').

Using the previous ISI sample and the best number of bins we construct the histogram and save it in file `u1hist1.txt` with:

    ./aspa_hist -n 46 -p --log < u1isi1.txt > u1hist1.txt

    Sample size: 3303
    Using a log transformation of the data.

The plot of the histogram with `gnuplot` is then done as follows:

    set grid
    unset key
    set title 'ISI density estimate neuron 1 spont. 1'
    set xlabel 'Inter spike interval (s)'
    set ylabel 'Estimated density (1/s)'
    plot [0:0.5] [] 'u1hist1.txt' using 1:3 with steps lc 'black'\
         linewidth 2

![img](fig/locust20010214_Spontaneous_1_tetB_u1_isi_hist.png)


<a id="orgd170863"></a>

## ISI histograms of the seven "good" neuron of dataset `locust20010214`

We can proceed and compute the ISI histograms of the 6 other good neurons of the experiment. We start by downloading the data with:

    prefix=https://raw.githubusercontent.com/christophe-pouzat/zenodo-locust-datasets-analysis
    prefix="$prefix"/master/Locust_Analysis_with_R/locust20010214/locust20010214_spike_trains/
    prefix="$prefix"locust20010214_Spontaneous_1_tetB_u
    for i in {2..7} ; do
        wget "$prefix"$i".txt"
    done

We now get the ISIs stored in files names `uXisi1.txt` where `X` takes successively values 2 to 7:

    file_name=locust20010214_Spontaneous_1_tetB_u
    for i in {2..7} ; do
        ./aspa_read_spike_train --inter_trial_interval=30 --trial_duration=29 < \
    			"$file_name"$i".txt" | \
    	./aspa_mst_isi > "u"$i"isi1.txt"
    done

We check if some ISI values are null:

    for i in {1..7} ; do
        in="u"$i"isi1.txt"
        sed -n '2,$p' < $in > dump
        echo "Unit" $i "has" $(grep "0$" < dump | wc -l) \
    	 "ISIs equal to zero among" $(($(wc -l < $in) - 1))
    done

    Unit 1 has 0 ISIs equal to zero among 3303
    Unit 2 has 0 ISIs equal to zero among 3574
    Unit 3 has 0 ISIs equal to zero among 1339
    Unit 4 has 0 ISIs equal to zero among 1890
    Unit 5 has 4 ISIs equal to zero among 4912
    Unit 6 has 0 ISIs equal to zero among 909
    Unit 7 has 1 ISIs equal to zero among 4155

Since these null ISIs are always making less than a thousands of the sample size, we remove them. We just have to be careful since the first line of the ISI files contains the number of events. It should therefore be changed if we remove lines:

    for i in {1..7} ; do
        # in contains input file name
        in="u"$i"isi1.txt"
        # remove first line of in and save result in dump
        sed -n '2,$p' < $in > dump
        # check how many lines with just 0 dump contains
        zero=$(grep "0$" < dump | wc -l)
        if [[ $zero != 0 ]] ; then
           # if there are lines with just zero remove them
           sed -i "/0$/d" $in
           # compute the new number of events in two steps
           original=$(sed -n '1p' < $in)
           new=$(($original - $zero))
           # update the number of events
           sed -i "1s/.*/$new/" $in
        fi
    done

We can then construct the ISI histograms (output not shown):

    for i in {1..7} ; do
        in="u"$i"isi1.txt"
        out="u"$i"_isi_dens1.txt"
        n=$(./aspa_hist_bw --best --log < $in)
        ./aspa_hist --n_bins=$n --log --prob < $in > $out
    done

We can now make a graph with all the estimated ISI densities as follows:

    set grid
    set title 'ISI density estimate for neurons 1 to 7 spont. 1'
    set xlabel 'Inter spike interval (s)'
    set ylabel 'Estimated density (1/s)'
    plot [0:0.5] [0:20] 'u1_isi_dens1.txt' using 1:3 with steps linewidth 2,\
         'u2_isi_dens1.txt' using 1:3 with steps linewidth 2,\
         'u3_isi_dens1.txt' using 1:3 with steps linewidth 2,\
         'u4_isi_dens1.txt' using 1:3 with steps linewidth 2,\
         'u5_isi_dens1.txt' using 1:3 with steps linewidth 2,\
         'u6_isi_dens1.txt' using 1:3 with steps linewidth 2,\
         'u7_isi_dens1.txt' using 1:3 with steps linewidth 2

![img](fig/locust20010214_Spontaneous_1_tetB_u1_7_isi_dens.png)

