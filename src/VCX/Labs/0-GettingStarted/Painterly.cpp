#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Labs/0-GettingStarted/Painterly.h"
#include "Labs/0-GettingStarted/PathTracing.h"
#include <iostream>
#include <functional>

namespace VCX::Labs::GettingStarted{

    //4. scene:radius, position, emission, color, material(0=diffuse,1=specular,2=refractive)
    Sphere spheres_painterly[]={
        Sphere(1e5, Vec(1e5+1, 40.8, 81.6), Vec(),Vec(0.75,0.25,0.25),0), //left wall
        Sphere(1e5, Vec(-1e5+99, 40.8, 81.6), Vec(),Vec(0.25,0.25,0.75),0), //right wall
        //Sphere(1e5, Vec(50, 40.8, 1e5), Vec(),Vec(0.75,0.75,0.75),0), //back wall
        Sphere(1e5, Vec(50, 40.8, 1e5), Vec(),Vec(0.25,0.25,0.25),0), //back wall
        Sphere(1e5, Vec(50, 40.8, -1e5+170), Vec(),Vec(),0), //front wall
        //Sphere(1e5, Vec(50, 1e5, 81.6), Vec(),Vec(0.75,0.75,0.75),0), //bottom wall
        Sphere(1e5, Vec(50, 1e5, 81.6), Vec(),Vec(0.25,0.25,0.25),0), //bottom wall
        Sphere(1e5, Vec(50, -1e5+81.6, 81.6), Vec(),Vec(0.75,0.75,0.75),0), //top wall
        Sphere(16.5, Vec(27, 16.5, 47), Vec(),Vec(1,1,1)*0.999,1), //mirror
        Sphere(16.5, Vec(73, 16.5, 78), Vec(),Vec(1,1,1)*0.999,2), //glass
        //Sphere(1.5, Vec(50, 81.6-16.5, 81.6), Vec(4,4,4)*100,Vec(),0), //light
        Sphere(600, Vec(50, 681.6-0.27, 81.6), Vec(12,12,12),Vec(),0), //light
    };
    int numSpheres_painterly = sizeof(spheres_painterly)/sizeof(Sphere);

    //painterly1: hemisphere
    Vec hemisphere(double u1, double u2) {
	    const double r=sqrt(1.0-u1*u1);
	    const double phi=2*M_PI*u2;
	    return Vec(cos(phi)*r,sin(phi)*r,u1);
    }

    //！！compute radiance estimate along ray (return radiance estimate)
    Vec radiance_painterly(const Ray &r, int depth, unsigned short *Xi, Halton& hal, Halton& hal2){
        //intersection
        double t; //最近交点距离
        int id = 0;//球的id
        if (!intersect(r, spheres_painterly,numSpheres_painterly,t,id)) return Vec() ;//miss:return black
        const Sphere &obj = spheres_painterly[id];//球的id
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
                return obj.e; //E:emission(?)
        }

    //0. diffuse reflection
    if(obj.refl == 0){
        //double r1 = 2 * M_PI * erand48(Xi);
        //double r2 = erand48(Xi), r2s = sqrt(r2); 
        //create orthonormal coordinate (w,u,v)
        Vec w = nl;
        Vec u = ((fabs(w.x)>0.1? Vec(0,1) : Vec(1))%w).norm();
        Vec v = w%u;
        //Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); //采样方向
            
        //painterly-3 
        hal.next();
        hal2.next();
        Vec d = (nl + hemisphere(hal.get(),hal2.get()));
        return obj.e + f.mult(radiance_painterly(Ray(x,d),depth,Xi,hal,hal2));
/*
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
                if (intersect(Ray(x,l),t,id) && id == i){
                    double omega = 2*M_PI*(1-cos_a_max); //计算立体角（solid angle）
                    e = e + f.mult(s.e*l.dot(nl)*omega)*M_1_PI;
                }
            }//loop over lights end
            return obj.e*E + e + f.mult(radiance(Ray(x,d),depth,Xi,hal, hal2,0));
*/
        }//diffuse end 

        //1. ideal specular(mirror)
        else if(obj.refl == 1){
            return obj.e + f.mult(radiance_painterly(Ray(x,r.d-n*2*n.dot(r.d)),depth, Xi, hal, hal2));
            //reflected ray: r.d-n*2*n.dot(r.d)
        }

        //2. dielectric (glass):reflective + refractive
        Ray reflRay(x, r.d-n*2*n.dot(r.d));
        bool into = n.dot(nl)>0; //ray from outside going in?
        //nc：起始介质折射率 nt：目标介质折射率（glass IOR=1.5） nnt：比值 ddn：入射角余弦
        double nc =1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;

        //判断是否发生全反射 cos2t=光线折射角的余弦平方值，cos2t<0全反射
        if((cos2t = 1-nnt*nnt*(1-ddn*ddn))<0){
            return obj.e + f.mult(radiance_painterly(reflRay, depth, Xi, hal, hal2));
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
            radiance_painterly(reflRay, depth, Xi,hal, hal2)*RP:radiance_painterly(Ray(x,tdir),depth,Xi,hal, hal2)*TP):
            radiance_painterly(reflRay, depth, Xi,hal, hal2)*Re+radiance_painterly(Ray(x,tdir),depth,Xi,hal, hal2)*Tr);

    }//Vec radiance end

    //loops over image pixels, creates image
    Vec* PathTracing_painterly(int w, int h, int samps, ProgressCallback progressCallback){

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

        #pragma omp parallel for schedule(dynamic, 1) private(r) //openmp
        //each loop should be run in its own thread

        //painterly
        // correlated Halton-sequence dimensions
	    Halton hal, hal2;
	    hal.number(0,2);
	    hal2.number(0,2);

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
                            r = r+radiance_painterly(Ray(cam.o + d*140, d.norm()),0,Xi,hal,hal2)*(1./samps);
                        }//tent filter end
                        //camera rays are pushed forward to start in interior？？
                        c[i] = c[i]+Vec(clamp(r.x),clamp(r.y),clamp(r.z))*0.25; //color accumulation with gamma correction
                    }
                }
            }
        }//loop over all pixels end

        //write out
        for (int i = 0; i < w * h; i++) {
/*          // 输出 c[i] 的值
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
}