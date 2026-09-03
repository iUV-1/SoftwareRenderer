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
        TexFile = new TGAImage();
        TexFile->read_tga_file(texPath.c_str());
        TexFile->flip_vertically();
        if(normalPath != "") {
            NormalFile = new TGAImage();
            NormalFile->read_tga_file(normalPath.c_str());
            NormalFile->flip_vertically();
            UseNormal = true;
        }
        if(specPath != "") {
            SpecularFile = new TGAImage();
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

    // move operator
    Material(Material&& other) noexcept {
        other.SpecularFile = SpecularFile;
        other.TexFile = TexFile;
        other.NormalFile = NormalFile;

        SpecularFile = nullptr;
        TexFile = nullptr;
        NormalFile = nullptr;
    }

    // copy operator
    // deep copy for safe memory
    Material(const Material& other) {
        if(other.SpecularFile) {
            SpecularFile = new TGAImage(*other.SpecularFile);
        }

        if(other.TexFile) {
            TexFile = new TGAImage(*other.TexFile);
        }

        if (other.NormalFile) {
            NormalFile = new TGAImage(*other.NormalFile);
        }
    }

    Material operator= (const Material &other) {
        if(other.SpecularFile) {
            SpecularFile = new TGAImage(*other.SpecularFile);
        }

        if(other.TexFile) {
            TexFile = new TGAImage(*other.TexFile);
        }

        if (other.NormalFile) {
            NormalFile = new TGAImage(*other.NormalFile);
        }
        return *this;
    }
    TGAImage *NormalFile;
    TGAImage *TexFile;
    TGAImage *SpecularFile;
    bool UseNormal = false;
    bool UseSpecular = false;
};

#endif //SOFTWARERENDERER_MATERIAL_HPP
