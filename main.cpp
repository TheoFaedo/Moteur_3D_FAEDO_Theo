#include "tgaimage.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

float norme(std::vector<float> n){
    return std::sqrt(std::pow(n[0], 2) + std::pow(n[1], 2) + std::pow(n[2], 2));
}

int recv_vertice_index(std::string &line){
    int taille = line.size();

    std::string tmp = "";

    for(int i = 0; i < taille; i++){
        if(line[i] != '/'){
            tmp += line[i];
        }else{
            break;
        }
    }
    return std::stoi(tmp)-1;
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
 void draw_filled_triangle(int x0, int y0, int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color){
    //Mise à l'origine du triangle
    std::vector<int> b = { x1-x0, y1-y0 };
    std::vector<int> c = { x2-x0, y2-y0 };
    
    float sb = 1.;
    float sc = 1.;
    if(b[0]!=0 || b[1]!=0){
        sb = 1/(std::sqrt(std::pow(b[0], 2) + std::pow(b[1], 2)));
    }

    if(c[0]!=0 || c[1]!=0){
        sc = 1/(std::sqrt(std::pow(c[0], 2) + std::pow(c[1], 2)));
    }

    std::vector<std::vector<int>> points = {  };

    float i = 0.;
    while(i<=1.){
        float j = 0.;
        float ji = j+i;
        while(j<=1. && ji<=1. && ji >= 0.){
            if(i+j<=1.){
                points.push_back({(int)std::floor(i*b[0]+j*c[0]), (int)std::floor(i*b[1]+j*c[1])});
            }
            j+=sc;
        }
        i+=sb;
    }

    for(int i = 0; i < points.size(); i++){
        image.set(points[i][0]+x0, points[i][1]+y0, color);
    }
}

int objParser(std::string path, std::vector<std::vector<int>> &sommets, std::vector<std::vector<int>> &faces, std::vector<std::vector<float>> &normales, int width = 2048, int height = 2048){
    std::ifstream file(path);

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string label;

            iss >> label;

            // Parsing
            if (label == "v") {
                float x, y, z;
                iss >> x >> y >> z;
                sommets.push_back({(int)(x*(width/2)+width/2), (int)(y*(height/2)+height/2), 0});
            } else if (label == "f"){
                std::string v1, v2, v3;
                iss >> v1 >> v2 >> v3;
                faces.push_back({recv_vertice_index(v1), recv_vertice_index(v2), recv_vertice_index(v3)});
            } else if (label == "vn"){
                float x, y, z;
                iss >> x >> y >> z;
                normales.push_back({x, y, z});
            } else {
                if(label != "vt") std::cerr << "Error parsing line: " << line << std::endl;
            }
        }
        file.close();
    } else {
        std::cerr << "Unable to open file." << std::endl;
        return 1;
    }
}

TGAColor getColor(std::vector<float> n1, std::vector<float> n2, std::vector<float> n3, TGAColor color){
    std::vector<float> n = {n1[0]+n2[0]+n3[0], n1[1]+n2[1]+n3[1], n1[2]+n2[2]+n3[2]};
    float no = norme(n);
    n = {n[0]/no, n[1]/no, n[2]/no};

    std::vector<float> lightdir = {0, 0, 1};
    //std::cout << n[0] << " " << n[1] << " " << n[2] << std::endl;
    float dot = n[0]*lightdir[0] + n[1]*lightdir[1] + n[2]*lightdir[2];
    dot = std::max(0.4f, dot);
    

    return {
        (std::uint8_t)(dot*color.bgra[0]),
        (std::uint8_t)(dot*color.bgra[1]),
        (std::uint8_t)(dot*color.bgra[2]),
        color.bgra[3]
    };
}


int main()
{
    constexpr int width = 2048;
    constexpr int height = 2048;
    TGAColor white = {{255, 255, 255, 255}, 4};
    TGAColor red = {{0, 0, 255, 255}, 4};
    
    TGAImage image(width, height, TGAImage::RGB);

    std::vector<std::vector<int>> sommets;
    std::vector<std::vector<int>> faces;
    std::vector<std::vector<float>> normales;

    objParser("../obj/african_head/african_head.obj", sommets, faces, normales, width, height);   

    for(int i = 0; i<faces.size(); i++){
        draw_filled_triangle(sommets[faces[i][0]][0], sommets[faces[i][0]][1], sommets[faces[i][1]][0], sommets[faces[i][1]][1], sommets[faces[i][2]][0], sommets[faces[i][2]][1], image, getColor(normales[faces[i][0]], normales[faces[i][1]], normales[faces[i][2]],  {{(std::uint8_t)std::rand(), (std::uint8_t)std::rand(), (std::uint8_t)std::rand(), 255}, 4}));
    }

    image.write_tga_file("test_out.tga", true, false);
    return 0;
}