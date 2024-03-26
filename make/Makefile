FC = caf
FCFLAGS = -O3

# ## ------------------------------------------------- Main Objective ---------------------------------------------------------

# TARGET = ../run/tsunami
TARGET = test
.PHONY: all clean
all: $(TARGET)

# ## ------------------------------------------------- PATHS and ObJS ---------------------------------------------------------
VPATH = ../src:../src/src_global:../src/src_io:../src/src_parallel:../src/src_functions:../src/src_parameters:../src/src_dynamic

vpath 

# mod_diff.o mod_field.o mod_parallel.o wave.o constantes.o

# OBJS = mod_io.o 

MAIN = main.f90

# ## ------------------------------------------------- Compiling rules ---------------------------------------------------------

$(TARGET): $(MAIN) $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

$(OBJS): %.o : %.f90
	$(FC) $(FCFLAGS) -c $<


# ## --------------------------------------------------- relies  ---------------------------------------------------------
# 阻止.mod到.o的隐式规则
%.o: %.mod

mod_field.o: mod_field.f90 mod_io.o mod_parallel.o
wave.o: wave.f90 mod_field.o

clean:
	-$(RM) tsunami *.o *.mod $(TARGET)