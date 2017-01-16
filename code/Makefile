P=programe_name
OBJECTS=
CFLAGS += `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11 
LDLIBS = `pkg-config --libs gsl ` 

$(P): $(OBJECTS)

all : libaspa.a

libaspa.a : aspa_single.o
	ar cr libaspa.a aspa_single.o

aspa_single_test_objects=aspa_single_test.o
aspa_single_test : $(aspa_single_test_objects) libaspa.a
	cc $(aspa_single_test_objects) libaspa.a $(LDLIBS) -o aspa_single_test

aspa_single_test.o : aspa.h

aspa_single_testB_objects=aspa_single_testB.o
aspa_single_testB : $(aspa_single_testB_objects) libaspa.a
	cc $(aspa_single_testB_objects) libaspa.a $(LDLIBS) -o aspa_single_testB

aspa_single_testB.o : aspa.h

aspa_single_testC_objects=aspa_single_testC.o
aspa_single_testC : $(aspa_single_testC_objects) libaspa.a
	cc $(aspa_single_testC_objects) libaspa.a $(LDLIBS) -o aspa_single_testC

aspa_single_testC.o : aspa.h

.PHONY : clean
clean :
	rm -f libaspa.a \
	aspa_single.o \
	$(aspa_single_test_objects) aspa_single_test \
	$(aspa_single_testB_objects) aspa_single_testB \
	$(aspa_single_testC_objects) aspa_single_testC
