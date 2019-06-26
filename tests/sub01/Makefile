TEST_SUB01_OBJS = sub.o test_sub01.o
FC = gfortran

all: test_sub01

test_sub01: sub.o test_sub01.o
	$(FC) -o test_sub01 $(TEST_SUB01_OBJS)

%.o:%.f90
	$(FC) -c $<

clean:
	@rm -f *~
	@rm -f *.o
	@rm -f *.mod
