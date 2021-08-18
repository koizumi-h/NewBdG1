
OBJ_DIR   = ./obj
BIN_DIR = .
SRC_DIR = ./src
MPI_OBJ_DIR   = ./mpiobj
#OUT_DIR = ./jexdata

# commands _/_/_/_/_/_/_/_/_
RM = rm


NR_MOD=$(OBJ_DIR)/nrtype.o $(OBJ_DIR)/nrutil.o $(OBJ_DIR)/nr.o
NR_OPT_OBJ= $(OBJ_DIR)/mnbrak.o $(OBJ_DIR)/brent.o $(OBJ_DIR)/linmin.o $(OBJ_DIR)/frprmn.o
OBJ = $(NR_MOD) $(OBJ_DIR)/math_mod.o
OBJ += $(OBJ_DIR)/utility_mod.o
OBJ += $(OBJ_DIR)/parameter_mod.o
OBJ +=  $(OBJ_DIR)/path_list_mod.o
OBJ +=  $(OBJ_DIR)/int_list_mod.o
OBJ +=  $(OBJ_DIR)/hop_iterator_mod.o
OBJ +=  $(OBJ_DIR)/site_list_mod.o 
OBJ +=  $(OBJ_DIR)/visual_mod.o
OBJ +=  $(OBJ_DIR)/base_mod.o 
OBJ +=  $(OBJ_DIR)/angular_variable_mod.o
OBJ +=  $(OBJ_DIR)/energy_mod.o
#OBJ +=  $(OBJ_DIR)/cnode_iterator_mod.o
#OBJ +=  $(OBJ_DIR)/cnode_creator_mod.o
OBJ +=  $(OBJ_DIR)/circle_list_mod.o
OBJ +=  $(OBJ_DIR)/fchi_optimizer_mod.o
#OBJ +=  $(OBJ_DIR)/fchi_optimizer2_mod.o
#OBJ +=  $(OBJ_DIR)/fchi_optimizer_ex_mod.o
#OBJ +=  $(OBJ_DIR)/fchi_optimizer_plusA_mod.o
OBJ +=  $(OBJ_DIR)/fchi_optimizer_Aem_mod.o
OBJ +=  $(OBJ_DIR)/io_mod.o
OBJ +=  $(OBJ_DIR)/conjugate_gradient_mod.o
OBJ +=  $(OBJ_DIR)/conjugate_gradient_rsoc_mod.o
#OBJ +=  $(OBJ_DIR)/conjugate_gradient_fchi_mod.o
OBJ +=  $(OBJ_DIR)/conjugate_gradient_spin_mod.o
OBJ +=  $(OBJ_DIR)/wave_function_mod.o
OBJ +=  $(OBJ_DIR)/car_parrinello_mod.o 
#OBJ +=  $(OBJ_DIR)/car_parrinello_wf_mod.o
#OBJ +=  $(OBJ_DIR)/basic_current_mod.o
#OBJ +=  $(OBJ_DIR)/single_wf_mod.o
#OBJ +=  $(OBJ_DIR)/single_wf_chi_mod.o
#OBJ +=  $(OBJ_DIR)/spinwave_mod.o
#OBJ +=  $(OBJ_DIR)/cranking_model_xi1.o
#OBJ +=  $(OBJ_DIR)/arpes_mod.o
#OBJ +=  $(OBJ_DIR)/car_parrinello_fchi_mod.o 
#OBJ +=  $(OBJ_DIR)/charge_constraint_lambda_calc_mod.o 
OBJ +=  $(OBJ_DIR)/orthogonal_lambda_calc_mod.o 
#OBJ +=  $(OBJ_DIR)/car_parrinello_beta_mod.o 
#OBJ +=  $(OBJ_DIR)/car_parrinello_gamma_mod.o 
#OBJ +=  $(OBJ_DIR)/car_parrinello_chi_mod.o 
OBJ +=  $(OBJ_DIR)/winding_number_mod.o 
####$(OBJ_DIR)/monte_carlo_zeta_mod.o
OBJ +=  $(OBJ_DIR)/measurement_mod.o 
OBJ +=  $(OBJ_DIR)/diagonalize_mod.o 
#OBJ +=  $(OBJ_DIR)/spinwave_mod.o
#OBJ +=  $(OBJ_DIR)/monte_carlo_zeta_mpi_mod.o
#OBJ +=  $(OBJ_DIR)/monte_carlo_zeta_mpi_mod2.o
#OBJ +=  $(OBJ_DIR)/monte_carlo_zeta_mpi_mod3.o

MPI_OBJ = $(subst $(OBJ_DIR), $(MPI_OBJ_DIR), $(OBJ))
MPI_OBJ += $(MPI_OBJ_DIR)/monte_carlo_2lay_mpi_mod.o


FC=ifort
#FC=ifort-17.0.2.163
MFC=mpiifort

#mOBJ=$(OBJ) $(OBJ_DIR)/main.o
#mOBJ=$(OBJ) $(OBJ_DIR)/main2.o
#mOBJ=$(OBJ) $(OBJ_DIR)/main3.o

mwf1=$(OBJ) $(OBJ_DIR)/make_wave_function1.o
mwf2=$(OBJ) $(OBJ_DIR)/make_wave_function2.o
#mwf3=$(OBJ) $(OBJ_DIR)/make_wave_function3.o
#je=$(OBJ) $(OBJ_DIR)/mwf2_jandE.o
#be=$(OBJ) $(OBJ_DIR)/mwf2_be.o
jb=$(OBJ) $(OBJ_DIR)/mwf2_jb.o
#mwfbeta=$(OBJ) $(OBJ_DIR)/make_wave_function_beta.o
#mwfgamma=$(OBJ) $(OBJ_DIR)/make_wave_function_gamma.o
#mwfbetajande=$(OBJ) $(OBJ_DIR)/mwf_beta_jandE.o
#mwfgammajande=$(OBJ) $(OBJ_DIR)/mwf_gamma_jandE.o
#deltacurrent=$(OBJ) $(OBJ_DIR)/delta_current.o

#test=$(OBJ) $(OBJ_DIR)/test.o
#diag_test=$(OBJ) $(OBJ_DIR)/diag_test.o
#diag_test4=$(OBJ) $(OBJ_DIR)/diag_test4.o
diag_sb=$(OBJ) $(OBJ_DIR)/diag_sb.o
diag_sb_sc=$(OBJ) $(OBJ_DIR)/diag_sb_sc.o
#diag_bulk=$(OBJ) $(OBJ_DIR)/diag_bulk.o
single_shot=$(OBJ) $(OBJ_DIR)/diag_single_shot.o

arpes=$(OBJ) $(OBJ_DIR)/arpes.o
stm=$(OBJ) $(OBJ_DIR)/stm.o
stm_p=$(OBJ) $(OBJ_DIR)/stm_point.o
fermi=$(OBJ) $(OBJ_DIR)/fermi.o


go=$(MPI_OBJ) $(MPI_OBJ_DIR)/main2.o

#FFLAGS=-Wall -fbounds-check -Wuninitialized
# -ffpe-trap=invalid,zero,overflow -fbacktrace
#FFLAGS= -O0 -g -llapack -lblas
#FFLAGS= -O3 -mkl -traceback -heap-arrays
#FFLAGS= -O3 -xHOST -ipo -mkl
#FFLAGS= -O3 -xHOST -ipo -mkl -qopenmp
#FFLAGS= -O3 -xHOST -ipo -mkl -parallel
#FFLAGS= -O0 -g -mkl -traceback -heap-arrays
#FFLAGS = -O0 -g -mkl -p -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback
#MFLAGS= -O3 -mkl -heap-arrays
#MFLAGS= -O0 -g -mkl -p -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -heap-arrays
#FFLAGS= -llapack -lblas

MFLAGS = -O3 -xHOST -mkl

FFLAGS = -O3 -xHOST -mkl
#FFLAGS = -O3 -xHOST -mkl -heap-arrays
#FFLAGS = -O0 -g -mkl -p -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -heap-arrays

ifeq ($(OPT),0) 
FFLAGS = -O0 -g -mkl -p -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback
#FFLAGS = -O0 -g -mkl -p -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -heap-arrays
MFLAGS = -O0 -g -mkl -p -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback
endif
ifeq ($(OPT),-1) 
#FFLAGS = -O3 -xHOST -ipo -mkl -qopt-report-phase=vec
FFLAGS = -O3 -xHOST -mkl -qopt-report-phase=vec
endif
ifeq ($(OPT),1) 
#FFLAGS = -O3 -xHOST -ipo -mkl
FFLAGS = -O3 -xHOST -mkl
#FFLAGS = -O3 -xHOST -mkl -heap-arrays
MFLAGS = -O3 -xHOST -mkl
endif
ifeq ($(OPT),2) 
#FFLAGS = -O3 -xHOST -ipo -mkl -qopenmp
FFLAGS = -O3 -xHOST -mkl -qopenmp
#FFLAGS = -O3 -xHOST -mkl -qopenmp -heap-arrays
#FFLAGS = -O3 -xHOST -mkl -qopenmp -pg
#FFLAGS = -O3 -xHOST -mkl -qopenmp -pg -qopt-report-phase=vec
endif
ifeq ($(OPT),3) 
FFLAGS = -O5 -xHOST -ipo -mkl -qopenmp
endif

#OBJECTS = $(wildcard $(OBJ_DIR)/*.o)
#MODULES = $(wildcard $(OBJ_DIR)/*.mod)

FMOD = -module $(OBJ_DIR)/

MFMOD = -module $(MPI_OBJ_DIR)/


.SUFFIXES:
.SUFFIXES:.o .f90
.PHONY:clean cdat all
.f90.o:
	$(FC) -c $< $(FFLAGS)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) -c $< $(FFLAGS) -o $@ $(FMOD)
$(MPI_OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(MFC) -c $< $(MFLAGS) -o $@ $(MFMOD)
#	$(MFC) -c $< $(MFLAGS) -o $@ $(FMOD)

#	$(MFC) -c $< $(MFLAGS)
#	$(FC) -c $< -CB -fpe1 -O0 -pg -stand -llpack -lblas -traceback -g
#	$(FC) -c $< -warn all -O0 -r8 -CB -fpe0 -stand -pg -openmp -parallel
#go:$(mOBJ)
##	$(MFC) -o $@.exe $(mOBJ) $(MFLAGS)
#	$(MFC) -o $(BIN_DIR)/$@.exe $(mOBJ) $(MFLAGS) $(FMOD)
#	$(FC) -o $@ $(OBJ) -CB -fpe1 -O0 -pg -stand -llapack -lblas -traceback -g
#	$(FC) -o $@ $(OBJ) -CB -fpe1 -O3 -pg -stand -llapack -lblas -parallel
#all:mwf1 mwf2 je
#all:mwf1
#all:mwf2 dsb dsb_sc
all:mwf1 mwf2 dsb dsb_sc single_shot arpes stm fermi go
mwf1:mdir $(mwf1)
	$(FC) -o $(BIN_DIR)/$@.exe $(mwf1) $(FFLAGS) $(FMOD)
mwf2:mdir $(mwf2)
	$(FC) -o $(BIN_DIR)/$@.exe $(mwf2) $(FFLAGS) $(FMOD)
#mwf3:mdir $(mwf3)
#	$(FC) -o $(BIN_DIR)/$@.exe $(mwf3) $(FFLAGS) $(FMOD)
#beta:mdir $(mwfbeta)
#	$(FC) -o $(BIN_DIR)/$@.exe $(mwfbeta) $(FFLAGS) $(FMOD)
#gamma:mdir $(mwfgamma)
#	$(FC) -o $(BIN_DIR)/$@.exe $(mwfgamma) $(FFLAGS) $(FMOD)
#je:$(je) mdir
#	$(FC) -o $(BIN_DIR)/$@.exe $(je) $(FFLAGS) $(FMOD)
#be:$(be) mdir
#	$(FC) -o $(BIN_DIR)/$@.exe $(be) $(FFLAGS) $(FMOD)
jb:$(jb) mdir
	$(FC) -o $(BIN_DIR)/$@.exe $(jb) $(FFLAGS) $(FMOD)
#jebeta:mdir $(mwfbetajande)
#	$(FC) -o $(BIN_DIR)/$@.exe $(mwfbetajande) $(FFLAGS) $(FMOD)
#jegamma:mdir $(mwfgammajande)
#	$(FC) -o $(BIN_DIR)/$@.exe $(mwfgammajande) $(FFLAGS) $(FMOD)
#deltacurrent:$(deltacurrent) mdir
#	$(FC) -o $(BIN_DIR)/$@.exe $(deltacurrent) $(FFLAGS) $(FMOD)
#test:mdir $(test)
#	$(FC) -o $(BIN_DIR)/$@.exe $(test) $(FFLAGS) $(FMOD)
#dtest:mdir $(diag_test)
#	$(FC) -o $(BIN_DIR)/$@.exe $(diag_test) $(FFLAGS) $(FMOD)
#dtest4:mdir $(diag_test4)
#	$(FC) -o $(BIN_DIR)/$@.exe $(diag_test4) $(FFLAGS) $(FMOD)
dsb:mdir $(diag_sb)
	$(FC) -o $(BIN_DIR)/$@.exe $(diag_sb) $(FFLAGS) $(FMOD)
dsb_sc:mdir $(diag_sb_sc)
	$(FC) -o $(BIN_DIR)/$@.exe $(diag_sb_sc) $(FFLAGS) $(FMOD)
single_shot:mdir $(single_shot)
	$(FC) -o $(BIN_DIR)/$@.exe $(single_shot) $(FFLAGS) $(FMOD)
arpes:mdir $(arpes)
	$(FC) -o $(BIN_DIR)/$@.exe $(arpes) $(FFLAGS) $(FMOD)
stm:mdir $(stm)
	$(FC) -o $(BIN_DIR)/$@.exe $(stm) $(FFLAGS) $(FMOD)
stm_p:mdir $(stm_p)
	$(FC) -o $(BIN_DIR)/$@.exe $(stm_p) $(FFLAGS) $(FMOD)
fermi:mdir $(fermi)
	$(FC) -o $(BIN_DIR)/$@.exe $(fermi) $(FFLAGS) $(FMOD)
#dbulk:mdir $(diag_bulk)
#	$(FC) -o $(BIN_DIR)/$@.exe $(diag_bulk) $(FFLAGS) $(FMOD)

go:mdir $(go)
	$(MFC) -o $@.exe $(go) $(MFLAGS) $(MFMOD)

clean:
	$(RM) -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod 2> /dev/null
	$(RM) -f $(MPI_OBJ_DIR)/*.o $(MPI_OBJ_DIR)/*.mod 2> /dev/null
#	rm -f *.o *.mod 2> /dev/null
cdat:
	rm -f *.fbin *.dat *.gp *.txt *.plt fort.* *.csv 2> /dev/null
cfig:
	rm -f *.eps *.ps *.jpg *.gif *.png 2> /dev/null

mdir:
	@if [ ! -d $(OBJ_DIR) ]; then \
		mkdir -p $(OBJ_DIR); \
		echo 'mkdir -p $(OBJ_DIR)'; \
  fi
	@if [ ! -d $(BIN_DIR) ]; then \
		mkdir -p $(BIN_DIR); \
		echo 'mkdir -p $(BIN_DIR)'; \
  fi
#	@if [ ! -d $(OUT_DIR) ]; then \
#		mkdir -p $(OUT_DIR); \
#	 echo 'mkdir -p $(OUT_DIR)'; \
#  fi
	@if [ ! -d $(MPI_OBJ_DIR) ]; then \
		mkdir -p $(MPI_OBJ_DIR); \
		echo 'mkdir -p $(MPI_OBJ_DIR)'; \
  fi
