# Moteur_3D_FAEDO_Theo

## Auteur
Théo FAEDO

## Présentation du programme

Ce programme est un petit moteur 3D permettant rendre un objet 3D texturés (à partir d'un fichier .obj et d'une texture tga).

Le moteur est conçu pour fonction via le CPU et n'est donc pas adapté pour du temps réel.

### Fonctionnalités

- lecture de fichiers .obj avec leur texture associée
- éclairage ambiant
- ombrage diffus
- affichage des triangles du modèle
- texture 
- utilisation d'un z-buffer
- caméra perspective
- choix de l'angle de vue

## Compiler et exécuter

Pour compiler ce programme, vous aurez besoin de cmake (version 3.28 de préférence).

Si vous n'avez pas le projet en local, clonez le projet via git :
```bash
>> git clone https://github.com/TheoFaedo/Moteur_3D_FAEDO_Theo.git
```

Rendez-vous à la racine du projet puis créer le repertoire *build* s'il n'éxiste pas déjà.
```bash
>> mkdir build && cd build
```

Une fois dans le repertoire build, compilez le projet dans celui-ci avec la commande suivante :
```bash
>> cmake .. && cmake --build .
```

Vous pouvez désormais exécuter le programme via la commande :
```bash
>> ./moteur3D <chemin_fichier_obj> <chemin_texture> <float_angle>
```

Un exemple de modèle et de texture est disponible dans le repertoire */obj*.
```bash
>> ./moteur3D "../obj/african_head/african_head.obj" "../obj/african_head/african_head_diffuse.tga" 0.5
```
Après exécution du programme, un fichier *build/out.tga* contient le résultat du rendu fait par le programme.

