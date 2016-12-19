# aspa
Analyse des séquences de potentiels d'action / Spike train analysis

## Analyse des séquences de potentiels d'action (de façon élémentaire)

L'idée de base est de mettre en œuvre la « philosophie » Unix/Linux -- telle qu'elle est par exemple exposée de façon spécifique dans l'article d'Arnold Robbins [_What's GNU_](http://www.linuxjournal.com/article/2762) et de façon plus générale dans [l'article wikipédia](https://fr.wikipedia.org/wiki/Philosophie_d%27Unix) -- pour l'analyse des séquences de potentiels d'action. Comme les dites séquences ne sont jamais trop gourmandes en mémoire, elles peuvent être stockées sous forme de fichier texte (ASCII) et le gros des opérations sur les séquences peuvent être conçues comme des « filtres » ; comprenez des programmes (en général mais pas forcément écris en `C`) qui lisent une entrée au format texte depuis « l'entrée standard » (`stdin`) et envoient leur résultat au format texte à « la sortie standard » (`stdout`). Pour les graphiques, [gnuplot](http://gnuplot.info/) va être utilisé.

## Spike train analysis (in an elementary way)

The idea here is to implement the Unix/Linux "philosophy"--as exposed for instance in the article of Arnold Robbins [_What's GNU_](http://www.linuxjournal.com/article/2762)--to the analysis of neuronal spike trains. Since spike trains make not too voluminous data, they can be stored as text files (ASCII) and most operations on them can be designed as "filters", that is programs (usually but not always written in `C`) that read their input in text format from the "standard input" (`stdin`) and send their result in text format to the "standard output" (`stdout`). For graphical displays, we are going to use [gnuplot](http://gnuplot.info/).
