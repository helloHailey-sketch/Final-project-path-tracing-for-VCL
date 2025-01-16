#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Labs/0-GettingStarted/PathTracing.h"
#include <iostream>
#include <functional>

namespace VCX::Labs::GettingStarted{
    //preliminaries
    //M_PI, M_1_PI=1/M_PI
    double erand48(unsigned short xsubi[3]){
        return (double)rand()/(double)RAND_MAX;
    }
/* 头文件
    //1. 设置Vec操作（用于坐标、向量、rgb）
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

    //2. ray structure
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
*/
    //4. scene:radius, position, emission, color, material(0=diffuse,1=specular,2=refractive)
    Sphere spheres[]={
        Sphere(1e5, Vec(1e5+1, 40.8, 81.6), Vec(),Vec(0.75,0.25,0.25),0), //left wall
        Sphere(1e5, Vec(-1e5+99, 40.8, 81.6), Vec(),Vec(0.25,0.25,0.75),0), //right wall
        Sphere(1e5, Vec(50, 40.8, 1e5), Vec(),Vec(0.75,0.75,0.75),0), //back wall
        Sphere(1e5, Vec(50, 40.8, -1e5+170), Vec(),Vec(),0), //front wall
        Sphere(1e5, Vec(50, 1e5, 81.6), Vec(),Vec(0.75,0.75,0.75),0), //bottom wall
        Sphere(1e5, Vec(50, -1e5+81.6, 81.6), Vec(),Vec(0.75,0.75,0.75),0), //top wall
        Sphere(16.5, Vec(27, 16.5, 47), Vec(),Vec(1,1,1)*0.999,1), //mirror
        Sphere(16.5, Vec(73, 16.5, 78), Vec(),Vec(1,1,1)*0.999,2), //glass
        Sphere(1.5, Vec(50, 81.6-16.5, 81.6), Vec(4,4,4)*100,Vec(),0), //light
    };
    int numSpheres = sizeof(spheres)/sizeof(Sphere);

/*见头文件
    //5. convert colors to displayable range [0,255]
    //inline double clamp(double x){return x<0 ? 0 : x>1 ? 1 : x;}
    //convert float to int; gamma correction of 2.2
    //inline int toInt(double x){ return int(pow(clamp(x), 1/2.2)*255+0.5);}

    //6. intersects ray with scene
    inline bool intersect(const Ray &r, double &t, int &id){
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
*/

    //bvh加速：struct AABB
    AABB::AABB() : min(Vec(1e20, 1e20, 1e20)), max(Vec(-1e20, -1e20, -1e20)) {}
    AABB::AABB(const Vec &min_, const Vec &max_) : min(min_), max(max_) {}
    // 包围盒与光线的相交测试
    bool AABB::intersect(const Ray &r, double &tmin, double &tmax) const {
        for (int i = 0; i < 3; i++) {
            double invD = 1.0 / r.d[i];//光线与平面x=min[i]相交的t=（min[i]-r.o[i])/r.d[i]
            double t0 = (min[i] - r.o[i]) * invD;//ray=r.o+t*r.d,进入包围盒的t
            double t1 = (max[i] - r.o[i]) * invD;//离开包围盒的t
            if (invD < 0.0) std::swap(t0, t1); //若光线方向为负，交换t0，t1
            tmin = t0 > tmin ? t0 : tmin;//光线进入包围盒
            tmax = t1 < tmax ? t1 : tmax;//光线离开包围盒
            if (tmax <= tmin) return false;
        }
        return true; //所有轴都相交：光线与包围盒相交
    }
    // 包围盒合并成一个大包围盒
    AABB AABB::merge(const AABB &a, const AABB &b) {
        Vec merge_min (std::min(a.min.x, b.min.x),std::min(a.min.y, b.min.y),std::min(a.min.z, b.min.z)) ;
        Vec merge_max (std::max(a.max.x, b.max.x),std::max(a.max.y, b.max.y),std::max(a.max.z, b.max.z)) ;
        return AABB(merge_min, merge_max);
    }
    //struct BVHNode见h
    //bvh构建:indices(球体索引数组),spheres(球体数组),start, end(球体范围)
    BVHNode *buildBVH(std::vector<int> &indices, const std::vector<Sphere> &spheres, int start, int end) {
        BVHNode *node = new BVHNode();

        // 计算当前节点的包围盒
        AABB box;
        for (int i = start; i < end; i++) {
            const Sphere &s = spheres[indices[i]];
            Vec sphereMin = s.p - Vec(s.rad, s.rad, s.rad);
            Vec sphereMax = s.p + Vec(s.rad, s.rad, s.rad);
            box = AABB::merge(box, AABB(sphereMin, sphereMax));
        }
        node->box = box;

        // 如果是叶节点：不再划分，返回node
        if (end - start <= 1) {
            node->start = start;
            node->end = end;
            return node;
        }

        // 划分空间，按坐标轴排序（x轴）
        int axis = 0;
        //依据中心点对indices排序，范围为start～end
        //！！indice才是真正记录了排序后的球，不是sphere下标！！
        std::sort(indices.begin() + start, indices.begin() + end, [&](int a, int b) {
            return spheres[a].p[axis] < spheres[b].p[axis];
        });

        // 递归构建左右子树
        int mid = (start + end) / 2;
        node->left = buildBVH(indices, spheres, start, mid);
        node->right = buildBVH(indices, spheres, mid, end);
        return node;
    }

    //调试node创建
    void printAABB(const AABB &box) {
    std::cout << "AABB: min=(" 
              << box.min.x << ", " << box.min.y << ", " << box.min.z << "), max=("
              << box.max.x << ", " << box.max.y << ", " << box.max.z << ")" 
              << std::endl;
}

    void printBVHNode(const BVHNode *node, int depth) {
    std::cout << std::string(depth * 2, ' ') // 用缩进表示深度
              << "Node at depth " << depth 
              << " | Start: " << node->start 
              << ", End: " << node->end 
              << std::endl;

    printAABB(node->box);
}
    void traverseBVH(const BVHNode *node, int depth = 0) {
    if (!node) return;

    // 打印当前节点信息
    printBVHNode(node, depth);

    // 递归遍历子节点
    if (node->left) traverseBVH(node->left, depth + 1);
    if (node->right) traverseBVH(node->right, depth + 1);
}

    //bvh相交测试:indices(球体索引数组),spheres(球体数组),start,end(球体范围)
    bool intersectBVH(const Ray &r, const BVHNode *node, const std::vector<Sphere> &spheres, const std::vector<int> &indices,double &t, int &id) {
        double tmin = 0, tmax = 1e20;
        if (!node->box.intersect(r, tmin, tmax)) return false; // 光线未与包围盒相交

        if (node->left == nullptr && node->right == nullptr) {
            // 没有子节点：自己是叶节点，测试光线与球体相交
            bool hit = false;
            for (int i = node->start; i < node->end; i++) {
                int sphereIdx = indices[i]; 
                double d = spheres[sphereIdx].intersect(r);
                if (d && d < t) {
                    t = d;
                    id = sphereIdx;
                    hit = true;
                }
            }
/*            std::cout << "Intersecting BVH, ray=(" << r.o.x << ", " << r.o.y << ", " << r.o.z 
          << ") -> (" << r.d.x << ", " << r.d.y << ", " << r.d.z << ")" << std::endl;
*/
            //if(id != 0) std::cout<<id;
            return hit;
        }
        // 检查左右子树
        bool hitLeft = node->left && intersectBVH(r, node->left, spheres, indices,t, id);
        bool hitRight = node->right && intersectBVH(r, node->right, spheres, indices, t, id);
        return hitLeft || hitRight;//光线与左右子树中的任意一个几何体相交
    }

/*
    //bvh相交接口
    bool intersect(const Ray &r, double &t, int &id, const BVHNode *bvhRoot, const std::vector<Sphere> &spheres) {
        return intersectBVH(r, bvhRoot, spheres, t, id);
    }
*/

    //！！compute radiance estimate along ray (return radiance estimate)
    Vec radiance(const Ray &r, int depth, unsigned short *Xi, const BVHNode *bvhRoot, const std::vector<Sphere> &spheres,const std::vector<int> &indices,int E){
        //intersection
        double t = 1e20; //最近交点距离
        int id = 0;//球的id
        //if (!intersect(r,spheres,numSpheres,t,id)) return Vec() ;//miss:return black
        if (!intersectBVH(r, bvhRoot, spheres,indices,t,id)) return Vec() ;
        const Sphere &obj = spheres[id];//球的id
        if (depth >10) return Vec();
        Vec x = r.o+r.d * t; //交点
        Vec n = (x-obj.p).norm(); //normal
        Vec nl = n.dot(r.d)<0?n:n*-1;//区分光线进入表面还是离开
        Vec f = obj.c; //颜色

        //Russian roulette
        double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y :f.z; //max reflectivity
        //depth5之后再轮盘赌
        if(++depth>5||!p) {
            if (erand48(Xi)<p) //如果随机数小于概率p，继续递归；否则，终止递归
                f = f * (1/p); //放大当前路径贡献
            else 
                return obj.e*E; //E:emission(?)
        }

        //0. diffuse reflection
        if(obj.refl == 0){
            double r1 = 2 * M_PI * erand48(Xi);
            double r2 = erand48(Xi), r2s = sqrt(r2); 
            //create orthonormal coordinate (w,u,v)
            Vec w = nl;
            Vec u = ((fabs(w.x)>0.1? Vec(0,1) : Vec(1))%w).norm();
            Vec v = w%u;
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); //采样方向


            //loop over lights
            Vec e;
            for(int i = 0; i<numSpheres; i++){
                const Sphere &s = spheres[i];
                if(s.e.x <= 0 && s.e.y <=0 && s.e.z <= 0) continue; //skip non-lights

                //create random direction towards sphere
                Vec sw = s.p -x, su=((fabs(sw.x)>0.1? Vec(0,1) : Vec(1))%sw).norm(), sv=sw%su; //create coordinate system for sampling: sw,su,sv
                double cos_a_max = sqrt(1-s.rad*s.rad/(x-s.p).dot(x-s.p)); //determine max angle
                //calculate sample direction based on randam numbers
                double eps1 = erand48(Xi), eps2 = erand48(Xi);
                double cos_a = 1-eps1 + eps1*cos_a_max;
                double sin_a = sqrt(1-cos_a*cos_a);
                double phi = 2 * M_PI * eps2;
                Vec l = su*cos(phi)*sin_a + sv * sin(phi)*sin_a + sw*cos_a;
                l.norm();

                //create shadow ray
                //if (intersect(Ray(x,l),spheres,numSpheres,t,id) && id == i){
                double shadowT = 1e20;
                if (intersectBVH(Ray(x,l),bvhRoot, spheres,indices, shadowT,id) && id == i){
                    double omega = 2*M_PI*(1-cos_a_max); //计算立体角（solid angle）
                    e = e + f.mult(s.e*l.dot(nl)*omega)*M_1_PI;
                }
            }//loop over lights end
            return obj.e*E + e + f.mult(radiance(Ray(x,d),depth,Xi,bvhRoot,spheres,indices,0));

        }//diffuse end 

        //1. ideal specular(mirror)
        else if(obj.refl == 1){
            return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth, Xi,bvhRoot,spheres,indices));
            //reflected ray: r.d-n*2*n.dot(r.d)
        }

        //2. dielectric (glass):reflective + refractive
        Ray reflRay(x, r.d-n*2*n.dot(r.d));
        bool into = n.dot(nl)>0; //ray from outside going in?
        //nc：起始介质折射率 nt：目标介质折射率（glass IOR=1.5） nnt：比值 ddn：入射角余弦
        double nc =1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;

        //判断是否发生全反射 cos2t=光线折射角的余弦平方值，cos2t<0全反射
        if((cos2t = 1-nnt*nnt*(1-ddn*ddn))<0){
            return obj.e + f.mult(radiance(reflRay, depth, Xi,bvhRoot,spheres,indices));
        }

        //otherwise反射或折射
        Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); //折射光线方向
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b),c=1-(into?-ddn:tdir.dot(n)); //R0:垂直入射时的菲涅耳反射系数，c=1-cos(theta)
        double Re = R0 + (1-R0)*c*c*c*c*c; //Schlick近似公式求菲涅耳反射系数
        double Tr = 1-Re;//折射比例
        //定义概率p，调整权重
        double P=0.25+0.5*Re, RP=Re/P, TP=Tr/(1-P);
        //depth>2则轮盘赌，erand48(Xi)<P折射，否则反射
        return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?
            radiance(reflRay, depth, Xi,bvhRoot,spheres,indices)*RP:radiance(Ray(x,tdir),depth,Xi,bvhRoot,spheres,indices)*TP):
            radiance(reflRay, depth, Xi,bvhRoot,spheres,indices)*Re+radiance(Ray(x,tdir),depth,Xi,bvhRoot,spheres,indices)*Tr);
    }//Vec radiance end

//进度条：声明回调函数类型
using ProgressCallback = std::function<void(float)>;

//loops over image pixels, creates image
Vec* PathTracing(int w, int h, int samps, ProgressCallback progressCallback){
    //image size
    //w=512, h=384;
    //sampling 
    //int samps = argc == 2 ? atoi(argv[1])/4 : 1;
    //samps = 1;
    //camera position & direction
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm());
    //horizantal(x) direction increment (uses implicit 0 for y,z);0.5135 defines field of view angle
    Vec cx = Vec(w*.5135/h); 
    Vec cy = (cx%cam.d).norm()*0.5135;//vertical(y) direction increment (note cross product)
    Vec r; //临时颜色
    Vec *c = new Vec[w*h]; //image result

    for (int i = 0; i < w * h; i++) {
        c[i] = Vec(0, 0, 0); // 将每个元素归零
    }

    //bvh结构初始化
    std::vector<Sphere> spheresVec(spheres, spheres + numSpheres);
    std::vector<int> indices(numSpheres);//索引
    for (int i = 0; i < numSpheres; i++) indices[i] = i;
    // 构建 BVH：indices(球体索引数组),spheres(球体数组),start, end(球体范围)
    BVHNode *bvhRoot = buildBVH(indices, spheresVec, 0, numSpheres);
    std::cout << "BVH Tree Structure:" << std::endl;
    traverseBVH(bvhRoot);
    #pragma omp parallel for schedule(dynamic, 1) private(r) //openmp
    //each loop should be run in its own thread

    //loop over all image pixels
    for(int y = 0; y < h; y++){
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%", samps*4,100.*y/(h-1)); //打印进度
        
        //进度条：更新进度，通过回调函数通知外部
        if(progressCallback){
            progressCallback(static_cast<float>(y)/h);
        }

        unsigned short Xi[3] = {0,0,static_cast<unsigned short>(y * y * y)}; //random number
        for(unsigned short x = 0; x<w; x++){
            //for each pixel do 2*2 subsample, samps samples per subsample
            for (int sy = 0, i = (h-y-1)*w+x; sy<2; sy++){ //2*2 subpixel rows
                for(int sx=0; sx<2; sx++, r = Vec()){ //2*2 subpixel columns
                    for (int s = 0; s < samps; s++){
                        //tent filter
                        double r1 = 2 * erand48(Xi), dx=r1<1? sqrt(r1)-1: 1-sqrt(2-r1);
                        double r2 = 2 * erand48(Xi), dy=r2<1? sqrt(r2)-1: 1-sqrt(2-r2);
                        Vec d = cx * (((sx+0.5+dx)/2 + x)/w-0.5) + cy * (((sy+0.5+dy)/2 + y)/h-0.5) + cam.d;
                        r = r+radiance(Ray(cam.o + d*140, d.norm()),0,Xi,bvhRoot,spheresVec,indices)*(1./samps);
                    }//tent filter end
                    //camera rays are pushed forward to start in interior？？
                    c[i] = c[i]+Vec(clamp(r.x),clamp(r.y),clamp(r.z))*0.25; //color accumulation with gamma correction
                }
            }
        }
    }//loop over all pixels end

    //write out
    for (int i = 0; i < w * h; i++) {
/*        // 输出 c[i] 的值
    std::cout << "c[" << i << "] = (" 
              << c[i].x << ", " 
              << c[i].y << ", " 
              << c[i].z << ")\n";
*/
    }

    //进度条：最终进度为100%
    if(progressCallback){
        progressCallback(1.0f);
    }

    return c;
}//path tracing() end

}//namespace end

