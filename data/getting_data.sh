# The Fatt and Katz data Table A.1.1 pp. 252-253 of Cox and Lewis
# The Statistical Analysis of Series of Events (1966) are available
# From Larry Wasserman web site (the are expressed in seconds).
wget http://www.stat.cmu.edu/~larry/all-of-nonpar/=data/nerve.dat

# Traffic data of Bartlett are available from Jim Lindsey's web
# site
wget http://www.commanster.eu/books/stoch2.tgz
tar -xvzf stoch2.tgz ch4/vehicles.dat
mv ch4/vehicles.dat .
rmdir ch4
rm stoch2.tgz


