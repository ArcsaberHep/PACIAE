#COMPILER and OPTIONS
FC=gfortran
LD=gfortran
FC_FLAGS=-O

#SOURCE and PROGRAM
OBJ_A=coales_22a.o hadcas_22a.o main_22a.o p22a.o parcas_22a.o sfm_22a.o
OBJ_B=coales_22b.o hadcas_22b.o main_22b.o p22b.o parcas_22b.o parini_22b.o sfm_22b.o
OBJ_C=coales_22c.o hadcas_22c.o main_22c.o p22c.o parcas_22c.o parini_22c.o sfm_22c.o
PRO_A=paciae_22a
PRO_B=paciae_22b
PRO_C=paciae_22c

all:
	@echo "usage: make run_a/run_b/run_c --- run program a/b/c"
	@echo "       make clean             --- clean object files"
run_a: $(PRO_A)
	time ./$(PRO_A)
run_b: $(PRO_B)
	time ./$(PRO_B)
run_c: $(PRO_C)
	time ./$(PRO_C)

$(PRO_A):$(OBJ_A)
	$(LD) -o $(PRO_A) $(OBJ_A)
	@echo ---------------
	@echo build successfully!
$(PRO_B):$(OBJ_B)
	$(LD) -o $(PRO_B) $(OBJ_B)
	@echo ---------------
	@echo build successfully!
$(PRO_C):$(OBJ_C)
	$(LD) -o $(PRO_C) $(OBJ_C)
	@echo ---------------
	@echo build successfully!
		

$(OBJ_A):%.o:%.f
	$(FC) $(FC_FLAGS) -c -o $@ $<
$(OBJ_B):%.o:%.f	
	$(FC) $(FC_FLAGS) -c -o $@ $<
$(OBJ_C):%.o:%.f
	$(FC) $(FC_FLAGS) -c -o $@ $<
	
clean:
	rm -f $(OBJ_A)
	rm -f $(OBJ_B)
	rm -f $(OBJ_C)
	rm -f $(PRO_A)
	rm -f $(PRO_B)
	rm -f $(PRO_C)
