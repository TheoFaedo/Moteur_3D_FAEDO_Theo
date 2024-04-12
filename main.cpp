#include "tgaimage.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#define C_SEUIL -0.01

float norme(std::vector<float> n){
    return std::sqrt(std::pow(n[0], 2) + std::pow(n[1], 2) + std::pow(n[2], 2));
}

// Fonction pour calculer le déterminant d'une matrice 2x2
float determinant2x2(float a, float b, float c, float d) {
    return a * d - b * c;
}

// Fonction pour calculer le déterminant d'une matrice 3x3
float determinant3x3(float matrix[9]) {
    return matrix[0] * determinant2x2(matrix[4], matrix[5], matrix[7], matrix[8]) -
           matrix[1] * determinant2x2(matrix[3], matrix[5], matrix[6], matrix[8]) +
           matrix[2] * determinant2x2(matrix[3], matrix[4], matrix[6], matrix[7]);
}

// Fonction pour calculer l'inverse d'une matrice 3x3
int inverse3x3(float matrix[9], float inverse[9]) {
    float det = determinant3x3(matrix);
    if (std::abs(det) < 1e-6) {
        return 1;
    }
    
    inverse[0] = determinant2x2(matrix[4], matrix[5], matrix[7], matrix[8]) / det;
    inverse[1] = -determinant2x2(matrix[1], matrix[2], matrix[7], matrix[8]) / det;
    inverse[2] = determinant2x2(matrix[1], matrix[2], matrix[4], matrix[5]) / det;
    inverse[3] = -determinant2x2(matrix[3], matrix[5], matrix[6], matrix[8]) / det;
    inverse[4] = determinant2x2(matrix[0], matrix[2], matrix[6], matrix[8]) / det;
    inverse[5] = -determinant2x2(matrix[0], matrix[2], matrix[3], matrix[5]) / det;
    inverse[6] = determinant2x2(matrix[3], matrix[4], matrix[6], matrix[7]) / det;
    inverse[7] = -determinant2x2(matrix[0], matrix[1], matrix[6], matrix[7]) / det;
    inverse[8] = determinant2x2(matrix[0], matrix[1], matrix[3], matrix[4]) / det;

    return 0;
}


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

float profondeur_z(float z0, float z1, float z2){
    return (z0+z1+z2)/3;
}

TGAColor getColor(std::vector<float> n, std::vector<int> t1, std::vector<int> t2, std::vector<int> t3, float alpha, float beta, float gamma, TGAImage &textureMap){

    std::vector<float> lightdir = {0, 0, 1};
    float dot = n[0]*lightdir[0] + n[1]*lightdir[1] + n[2]*lightdir[2];
    dot = std::max(0.3f, dot);

    TGAColor color = textureMap.get(t1[0]*alpha + t2[0]*beta + t3[0]*gamma, t1[1]*alpha + t2[1]*beta + t3[1]*gamma);
    //TGAColor color = {{255, 255, 255, 255}, 4};

    return {
        (u_int8_t)floor(color[0] * dot),
        (u_int8_t)floor(color[1] * dot),
        (u_int8_t)floor(color[2] * dot),
        (u_int8_t)floor(color[3] * dot),
    };
}

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

void draw_triangle(int x0, int y0, int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    draw_line(x0, y0, x1, y1, image, color);
    draw_line(x1, y1, x2, y2, image, color);
    draw_line(x2, y2, x0, y0, image, color);
}

/**
 * @brief Draw a filled triangle
 */ 
void draw_filled_triangle(int x0, int y0, int x1, int y1, int x2, int y2, TGAImage &image,
std::vector<int> t1, std::vector<int> t2, std::vector<int> t3, TGAImage &textureMap,
std::vector<float> n1, std::vector<float> n2, std::vector<float> n3, int* zbuffer, 
float z0, float z1, float z2){

    if(x0 > image.width() || y0 > image.height() || x1 > image.width() || y1 > image.height() || x2 > image.width() || y2 > image.height()
    || x0 < 0 || y0 < 0 || x1 < 0 || y1 < 0 || x2 < 0 || y2 < 0)
    return;

    std::vector<int> coinGauche = {std::min(x0, std::min(x1, x2)), std::min(y0, std::min(y1, y2))};
    std::vector<int> coinDroit = {std::max(x0, std::max(x1, x2)), std::max(y0, std::max(y1, y2))};

    float mat[9] = { (float)x0, (float)x1, (float)x2, (float)y0, (float)y1, (float)y2, 1., 1., 1. };
    float matInverse[9];

    if(inverse3x3(mat, matInverse)) return;

    #pragma omp parallel for
    for(int y = coinGauche[1]; y < coinDroit[1]; y++){
        for(int x = coinGauche[0]; x <= coinDroit[0]; x++){
            float alpha = matInverse[0] * x + matInverse[1] * y + matInverse[2];
            float beta = matInverse[3] * x + matInverse[4] * y + matInverse[5];
            float gamma = matInverse[6] * x + matInverse[7] * y + matInverse[8];
            if((alpha >= C_SEUIL && beta >= C_SEUIL && gamma >= C_SEUIL) || (alpha <= C_SEUIL && beta <= C_SEUIL && gamma <= C_SEUIL) && (alpha != C_SEUIL || beta != C_SEUIL || gamma != C_SEUIL)){
                float z = alpha * z0 + beta * z1 + gamma * z2;
                if(zbuffer[x+y*image.width()] < z){
                    zbuffer[x+y*image.width()] = z;
                    std::vector<float> n = {n1[0]*alpha+n2[0]*beta+n3[0]*gamma, n1[1]*alpha+n2[1]*beta+n3[1]*gamma, n1[2]*alpha+n2[2]*beta+n3[2]*gamma};
                    float no = norme(n);
                    n = {n[0]/no, n[1]/no, n[2]/no};
                    
                    //if(n[2] < 0) return;
                    image.set(x, y, getColor(n, t1, t2, t3, alpha, beta, gamma, textureMap));  
                }  
            }              
        }
    }
}

int objParser(std::string path, std::vector<std::vector<float>> &sommets, std::vector<std::vector<std::vector<int>>> &faces, std::vector<std::vector<float>> &normales, std::vector<std::vector<int>> &textures, TGAImage &textureMap, int width = 2048, int height = 2048){
    std::ifstream file(path);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string label;

            iss >> label;

            double angle = 4;

            // Parsing
            if (label == "v") {
                float x, y, z;
                iss >> x >> y >> z;
                //sommets.push_back({(int)(x*(width/2)+width/2), (int)(y*(height/2)+height/2), (int)(z*height)});
                sommets.push_back({(float)(x*cos(angle) - z*sin(angle)), y, (float)(x*sin(angle) + z*cos(angle))});
            } else if (label == "f"){
                std::string v1, v2, v3;
                iss >> v1 >> v2 >> v3;
                faces.push_back({{recv_vertice_index(v1), recv_vertice_index(v1, 1)}, {recv_vertice_index(v2), recv_vertice_index(v2, 1)}, {recv_vertice_index(v3), recv_vertice_index(v3, 1)}});
            } else if (label == "vn"){
                float x, y, z;
                iss >> x >> y >> z;
                normales.push_back({(float)(x*cos(angle) - z*sin(angle)), y, (float)(x*sin(angle) + z*cos(angle))});
            } else if (label == "vt"){
                float x, y, z;
                iss >> x >> y >> z;
                textures.push_back({(int)floor(x*textureMap.width()), textureMap.height()-((int)floor(y*textureMap.height())), 0});
            } else {
                if(label != "vt") std::cerr << "Error parsing line: " << line << std::endl;
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file." << std::endl;
        return 1;
    }

    return 0;
}

std::vector<int> to_perspective(std::vector<float> sommet, double tanAby2, int width, int height){

    float newX = sommet[0]/(((sommet[2]-3))*tanAby2);
    float newY = -sommet[1]/(((sommet[2]-3))*tanAby2);

    return {
        (int)(newX*(width/2)+width/2),
        (int)(newY*(height/2)+height/2),
        (int)(sommet[2]*height)
    };
}


int main()
{
    int model = 0;

    constexpr int width = 2048;
    constexpr int height = 2048;
    TGAColor white = {{255, 255, 255, 255}, 4};
    TGAColor red = {{0, 0, 255, 255}, 4};

    std::string obj_path;
    std::string texture_path;

    if(model == 0){
        obj_path = "../obj/african_head/african_head.obj";
        texture_path = "../obj/african_head/african_head_diffuse.tga";
    }else if(model == 1){
        obj_path = "../obj/boggie/body.obj";
        texture_path = "../obj/boggie/body_diffuse.tga";
    }else{
        obj_path = "../obj/diablo3_pose/diablo3_pose.obj";
        texture_path = "../obj/diablo3_pose/diablo3_pose_diffuse.tga";
    }
    
    TGAImage image(width, height, TGAImage::RGB);

    int* zbuffer = (int*)malloc(width*height*sizeof(int));
    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = -100000;
    }

    //Loading texture map
    TGAImage textureMap = TGAImage();
    textureMap.read_tga_file(texture_path);

    std::vector<std::vector<float>> sommetsRelative;
    std::vector<std::vector<int>> sommets;
    std::vector<std::vector<std::vector<int>>> faces;
    std::vector<std::vector<float>> normales;
    std::vector<std::vector<int>> textures;

    objParser(obj_path, sommetsRelative, faces, normales, textures, textureMap, width, height);

    double fov = M_PI/4.;
    double tanAby2 = tan(fov/2);
    std::cout << tanAby2 << std::endl;

    for(int i = 0; i < sommetsRelative.size(); i++){
        sommets.push_back(to_perspective(sommetsRelative[i], tanAby2, width, height));
        //std::cout << sommets[i][0] << " " << sommets[i][1] << " " << sommets[i][2] << std::endl;
    } 

    std::cout << faces.size() << std::endl;
    for(int i = 0; i<faces.size(); i++){
        //std::cout << i << std::endl;
        //std::cout << sommets[faces[i][0][0]][0] << " " << sommets[faces[i][0][0]][1] << " " << sommets[faces[i][0][0]][2] << std::endl;
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

    image.write_tga_file("test_out.tga", true, false);
    return 0;
}