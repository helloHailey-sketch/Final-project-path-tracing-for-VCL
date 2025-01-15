#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <functional>


namespace VCX::Labs::GettingStarted {

    // 1. Vec 类（用于坐标、向量、rgb）
    struct Vec{
        double x,y,z;  

        Vec(double x_=0, double y_=0, double z_=0){x=x_;y=y_;z=z_;}
        Vec operator+(const Vec &b) const {return Vec(x+b.x, y+b.y,z+b.z);}
        Vec operator-(const Vec &b) const {return Vec(x-b.x, y-b.y,z-b.z);}
        Vec operator*(double b) const {return Vec(x*b, y*b,z*b);} //&b引用const出错
        Vec mult(const Vec &b) const {return Vec(x*b.x, y*b.y, z*b.z);}
        Vec& norm() {return *this = *this *(1/sqrt(x*x+y*y+z*z));}
        double dot(const Vec &b)const {return x*b.x+y*b.y+z*b.z;}
        Vec operator%(Vec &b) {return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}//cross
    };

    // 2. Ray 类
    struct Ray{
        Vec o, d;
        //constructor（direction should always be normalized)
        Ray(Vec o_, Vec d_) : o(o_), d(d_){}
    };

/*
    // 3. Sphere 类
    struct Sphere {
        double rad; // 半径radius
        Vec p, e, c; // position（sphere center）,emission,color
        int refl; // reflection type: 0=diffuse,1=specular,2=refractive
        //constructor
        Sphere(double rad_, Vec p_, Vec e_, Vec c_, int refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
        
        //intersect(returns smallest distance, 0=nohit)
        double intersect(const Ray &r) const; //见cpp
    };
*/
    //5. convert colors to displayable range [0,255]
    inline double clamp(double x){return x<0 ? 0 : x>1 ? 1 : x;}
    //convert float to int; gamma correction of 2.2
    //inline int toInt(double x){ return int(pow(clamp(x), 1/2.2)*255+0.5);}
    inline double correct(double x){ return pow(clamp(x), 1/2.2);}
    
    
    //6. intersects ray with scene
    inline bool intersect(const Ray& r, double& t, int& id);

    // 4. 核心路径追踪函数
    Vec radiance(const Ray &r, int depth, unsigned short *Xi, int E = 1);

    //进度条：渲染进度回调函数类型
    using ProgressCallback = std::function<void(float)>;

    // 5. 路径追踪入口
    Vec* PathTracing(int w, int h, int samps,ProgressCallback progressCallback);
}
