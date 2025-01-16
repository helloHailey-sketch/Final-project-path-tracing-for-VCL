#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <functional>
#include <Labs/0-GettingStarted/PathTracing.h>

namespace VCX::Labs::GettingStarted {
    //painterly
    Vec hemisphere(double u1, double u2);
    //为什么声明显示无法重载orz

    //painterly2: hal
    class Halton {
	    double value, inv_base;
	    public:
	    void number(int i,int base) {
	    	double f = inv_base = 1.0/base;
	    	value = 0.0;
	    	while(i>0) {
		    	value += f*(double)(i%base);
		    	i /= base;
		    	f *= inv_base;
		    }
	    }
	    void next() {
		    double r = 1.0 - value - 0.0000001;
		    if(inv_base<r) value += inv_base;
		    else {
			    double h = inv_base, hh;
			    do {hh = h; h *= inv_base;} while(h >=r);
			    value += hh + h - 1.0;
		    }
	    }
	    double get() { return value; }
    };

    // 4. 核心路径追踪函数
    Vec radiance_painterly(const Ray &r, int depth, unsigned short *Xi, Halton& hal, Halton& hal2);

    //进度条：渲染进度回调函数类型
    using ProgressCallback = std::function<void(float)>;

    // 5. 路径追踪入口
    Vec* PathTracing_painterly(int w, int h, int samps,ProgressCallback progressCallback);
}