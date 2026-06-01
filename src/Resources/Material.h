//
// Created by iUV on 5/27/2026.
//

#ifndef SOFTWARERENDERER_MATERIAL_HPP
#define SOFTWARERENDERER_MATERIAL_HPP
#include "string"
#include "tgaimage.h"

class Material {
public:
    //Material(const Material *mat);
    Material() { }
    Material(std::string texPath, std::string normalPath = "", std::string specPath = "")
    {
        TexFile->read_tga_file(texPath.c_str());
        TexFile->flip_vertically();
        if(normalPath != "") {
            NormalFile->read_tga_file(normalPath.c_str());
            NormalFile->flip_vertically();
            UseNormal = true;
        }
        if(specPath != "") {
            SpecularFile->read_tga_file(specPath.c_str());
            SpecularFile->flip_vertically();
            UseSpecular = true;
        }
    }
    ~Material() {
        delete NormalFile;
        delete TexFile;
        delete SpecularFile;
    }
    TGAImage *NormalFile = new TGAImage();
    TGAImage *TexFile = new TGAImage();
    TGAImage *SpecularFile = new TGAImage();
    bool UseNormal = false;
    bool UseSpecular = false;
};



#endif //SOFTWARERENDERER_MATERIAL_HPP
