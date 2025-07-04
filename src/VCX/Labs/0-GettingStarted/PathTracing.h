#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <functional>
#include <algorithm>


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

        // aabb用：添加索引访问操作符
        double operator[](int i) const {
            if (i == 0) return x;
            if (i == 1) return y;
            if (i == 2) return z;
            throw std::out_of_range("Index out of range for Vec");
        }
    };

    // 2. Ray 类
    struct Ray{
        Vec o, d;
        //constructor（direction should always be normalized)
        Ray(Vec o_, Vec d_) : o(o_), d(d_){}
    };


    //3. sphere　
    struct Sphere{
        double rad; //radius
        Vec p, e, c; //position（sphere center）,emission,color
        int refl; //reflection type: 0=diffuse,1=specular,2=refractive

        //constructor
        Sphere(double rad_, Vec p_, Vec e_, Vec c_, int refl_): 
            rad(rad_),p(p_),e(e_),c(c_),refl(refl_){}

        //intersect(returns smallest distance, 0=nohit)
        double intersect(const Ray &r) const{
            //solve t^2*d.d+2*t*(o-p).d+(o-p).(o-p)-R^2=0
            Vec op = p - r.o;
            double t, eps = 1e-4;
            double b = op.dot(r.d);
            double det = b*b-op.dot(op)+rad*rad;
            if(det<0) return 0;
            else det = sqrt(det);
            return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
        }
    };
/*
    // 4. 声明 spheres 和 numSpheres
    extern Sphere spheres[];
    extern int numSpheres;
*/

    //bvh
    struct AABB {
        Vec min, max;
        AABB();
        AABB(const Vec &min_, const Vec &max_);
        bool intersect(const Ray &r, double &tmin, double &tmax) const;
        static AABB merge(const AABB &a, const AABB &b);
    };

    struct BVHNode {
        AABB box;          // 节点的包围盒
        int start, end;    // 球体的范围（叶节点）
        BVHNode *left, *right; // 子节点
        BVHNode() : start(0), end(0), left(nullptr), right(nullptr) {}
    };
    //bvh构建:indices(球体索引数组),spheres(球体数组),start, end(球体范围)
    BVHNode *buildBVH(std::vector<int> &indices, const std::vector<Sphere> &spheres, int start, int end);
    //光线与bvh相交
    bool intersectBVH(const Ray &r, const BVHNode *node, const std::vector<Sphere> &spheres, const std::vector<int> &indices,double &t, int &id);
    //使用bvh进行场景光线相交测试
    //bool intersect(const Ray &r, double &t, int &id, const BVHNode *bvhRoot, const std::vector<Sphere> &spheres);

    //5. convert colors to displayable range [0,255]
    inline double clamp(double x){return x<0 ? 0 : x>1 ? 1 : x;}
    //convert float to int; gamma correction of 2.2
    //inline int toInt(double x){ return int(pow(clamp(x), 1/2.2)*255+0.5);}
    inline double correct(double x){ return pow(clamp(x), 1/2.2);}
    
    
    //6. intersects ray with scene
    inline bool intersect(const Ray &r, const Sphere *spheres, int numSpheres, double &t, int &id){
        double d;
        double inf = t = 1e20;
        //double n = sizeof(spheres)/sizeof(Sphere);

        for(int i = numSpheres; i--;){
            if((d=spheres[i].intersect(r))&&d<t){
                t=d; //最近交点距离
                id=i; //球的id
            }
        }
        return t<inf; //true or false
    }


    // 4. 核心路径追踪函数
    Vec radiance(const Ray &r, int depth, unsigned short *Xi, const BVHNode *bvhRoot, const std::vector<Sphere> &spheres,const std::vector<int> &indices, int E = 1);

    //进度条：渲染进度回调函数类型
    using ProgressCallback = std::function<void(float)>;

    // 5. 路径追踪入口
    Vec* PathTracing(int w, int h, int samps,ProgressCallback progressCallback);
}
