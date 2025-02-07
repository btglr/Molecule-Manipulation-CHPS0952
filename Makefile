CC=gfortran
FLAGS=-fcheck=all -Wunused -pedantic -fimplicit-none

EXEC = ReadPDB RMSD
all: $(EXEC)

RMSD: transformation_obj.o vdw_obj.o vdw_manager.o atom_obj.o molecule_obj.o xyz_functions.o RMSD.o
	-@echo "Linking    $(@)"
	-@$(CC) $(FLAGS) -o $@.out $+

ReadPDB: transformation_obj.o vdw_obj.o vdw_manager.o atom_obj.o molecule_obj.o pdb_functions.o xyz_functions.o utilities.o ReadPDB.o
	-@echo "Linking    $(@)"
	-@$(CC) $(FLAGS) -o $@.out $+

%.o: %.f90
	-@echo ""
	-@echo "Generating $@"
	-@$(CC) $(FLAGS) -c $<

help:
	@echo "(C) JC.Boisson"
	@echo "Sous-commandes :"
	@echo "NomFichier (sans son extension)       : compile et fait l'édition de lien du fichier \"NomFichier.f90\" correspondant pour générer un exécutable"

	@echo "make cleanSource                      : supprime les fichiers parasites (*~, *.old,#*,*.bak)"
	@echo "make clean                            : supprime *tous* les fichiers reproductibles ici les .o, les .mod  et aussi les fichiers parasites"
	@echo "make clean_all                        : supprime *tous* les fichiers reproductibles, les fichiers parasites et aussi les exécutables"


###------------------------------
### Cleaning
###------------------------------------------------------------

clean:
	-@rm -rf *.o *.mod

clean_all: clean cleanSource
	-@rm -rf $(EXEC)

cleanSource:
	-@find . \( -name "*~" -o -name "*.old" -o -name "#*" -o -name "*.out" \) -print -exec rm \{\} \;


.PHONY:  $(EXEC) clean clean_all cleanSource
