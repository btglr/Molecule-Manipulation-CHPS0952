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

=== Utilisation du menu ===

Lors de l'exécution de ReadPDB, un menu interactif est présenté.
Celui-ci permet de choisir les transformations à appliquer ainsi que de vérifier la topologie de la molécule.
Pour chaque transformation, le choix est donné entre une transformation à valeur aléatoire ou à valeur choisie.
Si la valeur est insérée par l'utilisateur, un certain format doit être respecté.

Pour effectuer une translation d'un vecteur (à valeurs réelles) donné, celui-ci doit être séparé par des virgules.
Pour effectuer une rotation, la valeur de l'angle (à valeur réelle) doit être au format réel (exemple : 15. ou 15.0).

=== Commentaires ===

Le commentaire des fichiers XYZ contient sur les quatre premiers caractères le nombre de transformations effectuées, puis chaque transformation dans l'ordre dans lequel elle a été effectuée.
Ces transformations sont relues en mémoire à partir du fichier XYZ afin de pouvoir continuer à modifier une molécule qui aurait préalablement été écrite par ce même programme.