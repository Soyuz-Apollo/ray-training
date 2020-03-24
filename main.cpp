#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"


struct Svet
{
	Vec3f coord;
	float intensity;
	Svet(const Vec3f &c, const float &i): coord(c) , intensity(i){}
};

struct Material
{

	Vec3f color;
	Vec3f reflectivity;
	float reflExp;
	Material(const Vec3f &c, const Vec3f &r , const float &rE):color(c) , reflectivity(r), reflExp(rE){}
	Material(){}
};

struct Sphera
{
	Vec3f centr;
	float R;
	Material material;

	Sphera(const Vec3f &c, const float &rad , const Material &m): centr(c) , R(rad), material(m) {}

	bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &hit) const 
	{
		Vec3f L = centr - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > R*R) return false;
        float thc = sqrtf(R*R - d2);
        hit       = tca - thc;
        float t1 = tca + thc;
        if (hit < 0) hit = t1;
        if (hit < 0) return false;
        return true;
	}
};

Vec3f reflect(const Vec3f &I , const Vec3f &N)
{
	return I - N*2.f*(I*N);
}

bool first_hit(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphera> &spheri , Vec3f &hit, Vec3f &N, Material &material)
{
	float sfr_dist = std::numeric_limits<float>::max();
	float dist_i;

	for (int i = 0; i < spheri.size(); ++i)
	{
		if(spheri[i].ray_intersect(orig , dir , dist_i) && dist_i < sfr_dist)
		{
			sfr_dist = dist_i;
			hit = orig + dir*dist_i;
			N = (hit - spheri[i].centr).normalize();
			material = spheri[i].material;
		}
	}
	return sfr_dist < 1000;
}


Vec3f ray_start(const Vec3f &orig, const Vec3f &dir, const std::vector<Svet> svets, const std::vector<Sphera> spheri, int d = 0)
{
	float sphere_dist = std::numeric_limits<float>::max();
	float svet_diff = 0;
	float summ_col = 0;
	float svet_dist;
	float summ_refl = 0;
	Vec3f hit, Norm, svet_dir;
	Material material, matblya;
	Vec3f shadow_orig, shadow_pt, shadow_norm;
	Vec3f ref_dir, ref_orig, ref_col;

	
	if(d > 4 || !first_hit(orig , dir, spheri, hit , Norm ,material))
	{
		return Vec3f(0.2,0.7,0.8);
	}
	ref_dir = reflect(dir , Norm);
	ref_orig = ref_dir * Norm < 0 ? hit - Norm*1e-3 : hit + Norm*1e-3;
	ref_col = ray_start(ref_orig , ref_dir , svets, spheri , d+1);

	for (int i = 0; i < svets.size(); ++i)
	{
		svet_dir = (svets[i].coord - hit).normalize();
		svet_dist = (svets[i].coord - hit).norm();
		shadow_orig = svet_dir * Norm < 0 ? hit - Norm*1e-3 : hit + Norm*1e-3;
		if(first_hit(shadow_orig, svet_dir, spheri , shadow_pt , shadow_norm , matblya) && (shadow_pt - shadow_orig).norm() < svet_dist)
			continue;
		
		summ_col += svets[i].intensity * std::max(0.f , svet_dir*Norm);
		summ_refl += powf(std::max(0.f, -reflect(-svet_dir , Norm)*dir) , material.reflExp)*svets[i].intensity;
	}
	return material.color * summ_col* material.reflectivity[0] + Vec3f(1. , 1., 1.)*summ_refl *material.reflectivity[1] + ref_col * material.reflectivity[2];
}

void rendering()
{
	int WIDTH = 1024;
	int HEIGHT = 768;
	int fov = M_PI/2.;
	float x,y;
	Vec3f dir;

	std::vector<Vec3f> our_main_picture(WIDTH * HEIGHT);
	std::vector<Sphera> shari;
	std::vector<Svet> svets;

	shari.push_back(Sphera(Vec3f(0, 0, -16), 1.5 ,Material(Vec3f(0.4,0.4,0.3) ,Vec3f(1.f ,1.f ,0.0) , 50)));
	shari.push_back(Sphera(Vec3f(-5, 2,  -20), 3 ,Material(Vec3f(0.7,0.2,0.8) ,Vec3f(0.9 ,0.1 , 0.0) , 30)));
	shari.push_back(Sphera(Vec3f(-5, 4,  -12.3), 1.5 ,Material(Vec3f(0.0, 10.0, 0.8) ,Vec3f(0.6, 10.0, 0) ,1425)));
	shari.push_back(Sphera(Vec3f(5, 0,  -11), 1.0 ,Material(Vec3f(0.5, 0.0, 2) ,Vec3f(1.0 ,1.0, 0.) ,200)));
	shari.push_back(Sphera(Vec3f(-3, 0,  -11), 1.0 ,Material(Vec3f(0.5, 0.0, 2) ,Vec3f(1.0 ,0, 0.0) ,200)));
	shari.push_back(Sphera(Vec3f(1, -4,  -12), 1.0 ,Material(Vec3f(0.5, 0.0, 2) ,Vec3f(1.0 ,1.0, 0) ,200)));


	svets.push_back(Svet(Vec3f(-20, 20,  20), 1.0));
	svets.push_back(Svet(Vec3f(20,20,0) , 1.0));

	for (size_t j = 0; j<HEIGHT; j++) 
	{
        for (size_t i = 0; i<WIDTH; i++) 
        {
            x =  (2*(i + 0.5)/(float)WIDTH  - 1)*tan(fov/2.)*WIDTH/(float)HEIGHT;
            y = -(2*(j + 0.5)/(float)HEIGHT - 1)*tan(fov/2.);
            dir = Vec3f(x, y, -1).normalize();
            our_main_picture[i+j*WIDTH] = ray_start(Vec3f(0,0,0), dir, svets, shari );
        
        }
    }

	std :: ofstream result;

	result.open("./out.ppm");
	result << "P6\n" << WIDTH << " " << HEIGHT << "\n255\n";

    for (size_t i = 0; i < WIDTH*HEIGHT; ++i) {
        for (size_t j = 0; j<3; j++) {
            result << (char)(255 * std::max(0.f, std::min(1.f, our_main_picture[i][j])));
        }
    }
        result.close();

}

int main()
{
	rendering();
	return 0;
}