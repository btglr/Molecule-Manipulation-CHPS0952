=== Compilation ===

make

=== Utilisation ===

Deux exécutables sont présents :

ReadPDB.out
RMSD.out

Tous les fichiers PDB sont mis à disposition dans le dossier PDB.
Les fichiers de sortie seront générés dans le dossier XYZ.

=== ReadPDB ===

./ReadPDB.out pdb_file [number_of_output_files]

=== RMSD ===

./RMSD.out xyz_file_1 xyz_file_2

=== Exemples ===

./ReadPDB.out PDB/7cmm.pdb 5
./ReadPDB.out PDB/7kah.pdb

./RMSD.out XYZ/7cmm_0001.xyz XYZ/7cmm_0003.xyz