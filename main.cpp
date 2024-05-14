#include "tgaimage.h"
#include "geometry.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#define C_SEUIL -0.01 // Permet de corriger les problèmes d'espacements des triangles

/**
 * @brief Permet de calculer la norme d'un vecteur
*/
float norme(std::vector<float> n){
    return std::sqrt(std::pow(n[0], 2) + std::pow(n[1], 2) + std::pow(n[2], 2));
}

/**
 * @brief Recupere l'index d'un sommet d'une composante d'un vecteur de face
 * @example recv_vertice_index("1/5/4") == 1 and recv_vertice_index("1/5/4", 1) == 5
*/
int recv_vertice_index(std::string &line, int index = 0){
    int taille = line.size();

    std::string tmp = "";

    int ind = -1;

    for(int i = 0; i < taille; i++){
        if(line[i] != '/'){
            tmp += line[i];
        }else{
            ind++;
            if(ind == index){
                break;
            }else{
                tmp="";
            }
        }
    }
    return std::stoi(tmp)-1;
}

/**
 * @brief Permet de recuperer la couleur d'un pixel à partir de ses coordonnées barycentriques dans le triangle,
 * des vecteurs textures du triangle, de la texture et de la normale de la face.
 * 
 * @param n Vecteur normal de la face
 * @param t1, t2, t3 Vecteurs textures
 * @param alpha, beta, gamma Coordonnées barycentriques
 * @param textureMap Image de texture
*/
TGAColor getColor(std::vector<float> n, std::vector<int> t1, std::vector<int> t2, std::vector<int> t3, float alpha, float beta, float gamma, TGAImage &textureMap){

    // Récupération de l'illumination de la face à partir de la norme
    std::vector<float> lightdir = {0, 0, 1};
    float dot = n[0]*lightdir[0] + n[1]*lightdir[1] + n[2]*lightdir[2];
    dot = std::max(0.3f, dot); // Lumière ambiante de 0.3

    // Application de la couleur à partir de la texture
    TGAColor color = textureMap.get(t1[0]*alpha + t2[0]*beta + t3[0]*gamma, t1[1]*alpha + t2[1]*beta + t3[1]*gamma);

    return {
        (u_int8_t)floor(color[0] * dot),
        (u_int8_t)floor(color[1] * dot),
        (u_int8_t)floor(color[2] * dot),
        (u_int8_t)floor(color[3] * dot),
    };
}

/**
 * @brief Permet de dessiner une ligne
 * @param x0, y0, x1, y1 : coordonnées des deux extremites
 * @param image : image dans laquelle on dessine
 * @param color : couleur de la ligne
*/
void draw_line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color){
    bool steep = false; 
    if (std::abs(x0-x1)<std::abs(y0-y1)) { 
        std::swap(x0, y0); 
        std::swap(x1, y1); 
        steep = true; 
    } 
    if (x0>x1) { 
        std::swap(x0, x1); 
        std::swap(y0, y1); 
    } 
    int dx = x1-x0; 
    int dy = y1-y0; 
    int derror2 = std::abs(dy)*2; 
    int error2 = 0; 
    int y = y0; 
    for (int x=x0; x<=x1; x++) { 
        if (steep) { 
            image.set(y, x, color); 
        } else { 
            image.set(x, y, color); 
        } 
        error2 += derror2; 
        if (error2 > dx) { 
            y += (y1>y0?1:-1); 
            error2 -= dx*2; 
        } 
    } 
}

/**
 * @brief Permet de dessiner un triangle
 * @param x0, y0, x1, y1, x2, y2 : coordonnées des trois extremites
 * @param image : image dans laquelle on dessine
 * @param color : couleur du triangle
*/
void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    draw_line(x0, y0, x1, y1, image, color);
    draw_line(x1, y1, x2, y2, image, color);
    draw_line(x2, y2, x0, y0, image, color);
}

/**
 * @brief Draw a filled triangle
 * @param x0, y0, x1, y1, x2, y2 : coordonnées des trois extremites
 * @param image : image dans laquelle on dessine
 * @param t1, t2, t3 : vecteurs textures
 * @param n1, n2, n3 : vecteurs normaux
 * @param zbuffer : zbuffer
 * @param z0, z1, z2 : profondeur des trois extremites
 * @param textureMap : texture
 */ 
void draw_filled_triangle(int x0, int y0, int x1, int y1, int x2, int y2, TGAImage &image,
std::vector<int> t1, std::vector<int> t2, std::vector<int> t3, TGAImage &textureMap,
std::vector<float> n1, std::vector<float> n2, std::vector<float> n3, int* zbuffer, 
float z0, float z1, float z2){

    // On ne dessine que si les points sont dans l'image
    if(x0 > image.width() || y0 > image.height() || x1 > image.width() || y1 > image.height() || x2 > image.width() || y2 > image.height()
    || x0 < 0 || y0 < 0 || x1 < 0 || y1 < 0 || x2 < 0 || y2 < 0)
    return;

    // On calcule les coordonnées des extremites
    std::vector<int> coinGauche = {std::min(x0, std::min(x1, x2)), std::min(y0, std::min(y1, y2))};
    std::vector<int> coinDroit = {std::max(x0, std::max(x1, x2)), std::max(y0, std::max(y1, y2))};

    // Calcul de la matrice inverse (pour le calcul des coéficients barycentriques)
    mat<3,3> ma = {
        vec<3>{(double)x0, (double)x1, (double)x2},
        vec<3>{(double)y0, (double)y1, (double)y2},
        vec<3>{1., 1., 1.}
    };
    mat<3,3> matInverse = ma.invert();

    #pragma omp parallel for
    for(int y = coinGauche[1]; y < coinDroit[1]; y++){
        for(int x = coinGauche[0]; x <= coinDroit[0]; x++){

            // Calcul des coefficients barycentriques
            float alpha = matInverse[0][0] * x + matInverse[0][1] * y + matInverse[0][2];
            float beta = matInverse[1][0] * x + matInverse[1][1] * y + matInverse[1][2];
            float gamma = matInverse[2][0] * x + matInverse[2][1] * y + matInverse[2][2];
            
            // Dessin des points dans triangle
            if((alpha >= C_SEUIL && beta >= C_SEUIL && gamma >= C_SEUIL) || (alpha <= C_SEUIL && beta <= C_SEUIL && gamma <= C_SEUIL) && (alpha != C_SEUIL || beta != C_SEUIL || gamma != C_SEUIL)){
                
                // Calcul de la profondeur du point
                float z = alpha * z0 + beta * z1 + gamma * z2;

                // On ne dessine que si le points est moins profond que le précédement (s'il y en a un)
                if(zbuffer[x+y*image.width()] < z){
                    zbuffer[x+y*image.width()] = z;

                    // Calcul de la normale de la face
                    std::vector<float> n = {n1[0]*alpha+n2[0]*beta+n3[0]*gamma, n1[1]*alpha+n2[1]*beta+n3[1]*gamma, n1[2]*alpha+n2[2]*beta+n3[2]*gamma};
                    float no = norme(n);
                    n = {n[0]/no, n[1]/no, n[2]/no}; // Normalisation
                    
                    if(n[2] < 0) return; // Si la norme est dans la direction de la profondeur alors on ne dessine pas (peu provoquer des problèmes dans certains cas)

                    // Application de la couleur à l'image  
                    image.set(x, y, getColor(n, t1, t2, t3, alpha, beta, gamma, textureMap));
                }  
            }              
        }
    }
}

/**
 * @brief Permet de parser un fichier .obj pour l'utiliser dans le programme
 * @param path Chemin du fichier .obj
 * @param sommets, faces, normales, textures : vecteurs qui contiennent les informations de l'objet
 * @param textureMap : image de texture
 * @param angle : angle de rotation
*/
int objParser(std::string path, 
std::vector<std::vector<float>> &sommets, std::vector<std::vector<std::vector<int>>> &faces, std::vector<std::vector<float>> &normales, 
std::vector<std::vector<int>> &textures, 
TGAImage &textureMap, double angle = 0.4){
    std::ifstream file(path);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string label;

            iss >> label; // lecture du label (v, vn, f, vt)

            // Parsing
            if (label == "v") { // sommets
                float x, y, z;
                iss >> x >> y >> z;
                sommets.push_back({(float)(x*cos(angle) - z*sin(angle)), y, (float)(x*sin(angle) + z*cos(angle))});

            } else if (label == "f"){ // faces
                std::string v1, v2, v3;
                iss >> v1 >> v2 >> v3;
                faces.push_back({{recv_vertice_index(v1), recv_vertice_index(v1, 1)}, {recv_vertice_index(v2), recv_vertice_index(v2, 1)}, {recv_vertice_index(v3), recv_vertice_index(v3, 1)}});

            } else if (label == "vn"){ // normales
                float x, y, z;
                iss >> x >> y >> z;
                normales.push_back({(float)(x*cos(angle) - z*sin(angle)), y, (float)(x*sin(angle) + z*cos(angle))});

            } else if (label == "vt"){ // textures
                float x, y, z;
                iss >> x >> y >> z;
                textures.push_back({(int)floor(x*textureMap.width()), textureMap.height()-((int)floor(y*textureMap.height())), 0});
            
            } else {
                if(label != "#" && label != "" && label != "g" && label != "s") std::cerr << "Error parsing line: " << line << std::endl;
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file." << std::endl;
        return 1;
    }

    return 0;
}

/**
 * @brief Permet de transformer par la perspective un point de coordonnées orthographiques
 * @param sommet : point de coordonnées orthographiques
 * @param tanAby2 : tangente de l'angle du champs de vision divisé par deux
 * @param width, height : largeur et hauteur de l'image
*/
std::vector<int> to_perspective(std::vector<float> sommet, double tanAby2, int width, int height){

    float cameraZ = 3;

    float newX = sommet[0]/(((sommet[2]-cameraZ))*tanAby2);
    float newY = -sommet[1]/(((sommet[2]-cameraZ))*tanAby2);

    return {
        (int)(newX*(width/2)+width/2),
        (int)(newY*(height/2)+height/2),
        (int)(sommet[2]*height)
    };
}


int main(int argc, char *argv[])
{

    float angle = 0.0;
    std::string obj_path;
    std::string texture_path;

    if(argc >= 3 && argc <= 4){
        obj_path = argv[1];
        texture_path = argv[2];
        if(argc == 4) angle = std::stof(argv[3]);
    } else {
        std::cerr << "Usage: " << argv[0] << " <obj_path> <texture_path> <angle>" << std::endl;
        return 1;
    }

    constexpr int width = 2048;
    constexpr int height = 2048;
    
    // Initialisation de l'image
    TGAImage image(width, height, TGAImage::RGB);

    // Initialisation du zbuffer
    int* zbuffer = (int*)malloc(width*height*sizeof(int));
    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = -100000;
    }

    // Loading texture map
    TGAImage textureMap = TGAImage();
    textureMap.read_tga_file(texture_path);

    // Création des listes de sommets, normales, textures et faces
    std::vector<std::vector<float>> sommetsRelative;
    std::vector<std::vector<int>> sommets;
    std::vector<std::vector<std::vector<int>>> faces;
    std::vector<std::vector<float>> normales;
    std::vector<std::vector<int>> textures; 
    
    // Parsing du fichier obj
    objParser(obj_path, sommetsRelative, faces, normales, textures, textureMap, angle);

    // Création de la valeur de FOV
    double fov = M_PI/4.;
    double tanAby2 = tan(fov/2);

    // Transformation en vecteurs perspective
    for(int i = 0; i < sommetsRelative.size(); i++){
        sommets.push_back(to_perspective(sommetsRelative[i], tanAby2, width, height));
    } 

    // Dessin des triangles sur l'image
    for(int i = 0; i<faces.size(); i++){
        draw_filled_triangle(
            sommets[faces[i][0][0]][0], sommets[faces[i][0][0]][1], 
            sommets[faces[i][1][0]][0], sommets[faces[i][1][0]][1], 
            sommets[faces[i][2][0]][0], sommets[faces[i][2][0]][1],
            image, 
            textures[faces[i][0][1]], textures[faces[i][1][1]], textures[faces[i][2][1]], textureMap,
            normales[faces[i][0][0]], normales[faces[i][1][0]], normales[faces[i][2][0]], 
            zbuffer,
            sommets[faces[i][0][0]][2], sommets[faces[i][1][0]][2], sommets[faces[i][2][0]][2]
        );
    }

    //Sauvegarde de l'image
    image.write_tga_file("out.tga", true, false);
    return 0;
}