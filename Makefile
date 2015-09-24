EXEC	:= libdata_trans_lib.a
VPATH	:= .
CXX		:= mpiicpc
FC		:= mpiifort
LD		:= mpiifort
AR		:= ar
ARFLAGS	:= rv
SRCS	:= $(wildcard $(addsuffix /*.cxx, $(VPATH)) $(addsuffix /*.F90, $(VPATH)))
OBJS	:= $(addsuffix .o, $(sort $(basename $(notdir $(SRCS)))))

.SUFFIXES:
.SUFFIXES: .cxx .F90 .o

all: $(EXEC)

$(EXEC): $(OBJS)
	$(AR) $(ARFLAGS) $(EXEC) $(OBJS)

.F90.o:
	$(FC) -c -I. $<

.cxx.o:
	$(CXX) -c -I. $<

clean:
	rm *.o $(EXEC) *.mod
