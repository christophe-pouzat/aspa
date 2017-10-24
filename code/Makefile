P=programe_name
OBJECTS=
CFLAGS += `pkg-config --cflags gsl` -g -Wall -O0 -std=gnu11 
LDLIBS = `pkg-config --libs gsl ` 

$(P): $(OBJECTS)

all : libaspa.a aspa_read_spike_train aspa_mst_fns aspa_mst_aggregate aspa_mst_plot\
aspa_mst_isi aspa_hist_bw

libaspa.a : aspa_single.o
	ar cr libaspa.a aspa_single.o

aspa_read_spike_train_objects=aspa_read_spike_train.o
aspa_read_spike_train : $(aspa_read_spike_train_objects) libaspa.a
	cc $(aspa_read_spike_train_objects) libaspa.a $(LDLIBS) -o aspa_read_spike_train

aspa_read_spike_train.o : aspa.h

aspa_mst_fns_objects=aspa_mst_fns.o
aspa_mst_fns : $(aspa_mst_fns_objects) libaspa.a
	cc $(aspa_mst_fns_objects) libaspa.a $(LDLIBS) -o aspa_mst_fns

aspa_mst_fns.o : aspa.h

aspa_mst_aggregate_objects=aspa_mst_aggregate.o
aspa_mst_aggregate : $(aspa_mst_aggregate_objects) libaspa.a
	cc $(aspa_mst_aggregate_objects) libaspa.a $(LDLIBS) -o aspa_mst_aggregate

aspa_mst_aggregate.o : aspa.h

aspa_mst_plot_objects=aspa_mst_plot.o
aspa_mst_plot : $(aspa_mst_plot_objects) libaspa.a
	cc $(aspa_mst_plot_objects) libaspa.a $(LDLIBS) -o aspa_mst_plot

aspa_mst_plot.o : aspa.h

aspa_mst_isi_objects=aspa_mst_isi.o
aspa_mst_isi : $(aspa_mst_isi_objects) libaspa.a
	cc $(aspa_mst_isi_objects) libaspa.a $(LDLIBS) -o aspa_mst_isi

aspa_mst_isi.o : aspa.h

aspa_hist_bw_objects=aspa_hist_bw.o
aspa_hist_bw : $(aspa_hist_bw_objects) libaspa.a
	cc $(aspa_hist_bw_objects) libaspa.a $(LDLIBS) -o aspa_hist_bw

aspa_hist_bw.o : aspa.h

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

aspa_single_testD_objects=aspa_single_testD.o
aspa_single_testD : $(aspa_single_testD_objects) libaspa.a
	cc $(aspa_single_testD_objects) libaspa.a $(LDLIBS) -o aspa_single_testD

aspa_single_testD.o : aspa.h

aspa_single_testE_objects=aspa_single_testE.o
aspa_single_testE : $(aspa_single_testE_objects) libaspa.a
	cc $(aspa_single_testE_objects) libaspa.a $(LDLIBS) -o aspa_single_testE

aspa_single_testE.o : aspa.h

.PHONY : clean
clean :
	rm -f libaspa.a \
	aspa_single.o \
	$(aspa_read_spike_train_objects) aspa_read_spike_train \
	$(aspa_mst_fns_objects) aspa_mst_fns \
	$(aspa_mst_aggregate_objects) aspa_mst_aggregate \
	$(aspa_mst_plot_objects) aspa_mst_plot \
	$(aspa_mst_isi_objects) aspa_mst_isi \
	$(aspa_hist_bw_objects) aspa_hist_bw \
	$(aspa_single_test_objects) aspa_single_test \
	$(aspa_single_testB_objects) aspa_single_testB \
	$(aspa_single_testC_objects) aspa_single_testC \
	$(aspa_single_testD_objects) aspa_single_testD \
	$(aspa_single_testE_objects) aspa_single_testE
