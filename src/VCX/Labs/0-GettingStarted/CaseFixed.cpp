#include <algorithm>
#include <array>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "Labs/0-GettingStarted/CaseFixed.h"
#include "Labs/Common/ImGuiHelper.h"
namespace VCX::Labs::GettingStarted{
    //preliminaries
    //1. 设置vec计算操作（用于坐标、向量、rgb）
    struct Vec{
      double x,y,z;  

    };

    //2. ray structure
    struct Ray{
        Vec ray_origin, ray_direction;
        //constructor（direction should always be normalized)
        Ray(Vec ray_origin_, ray_direction_) : ray_origin(ray_origin_), ray_direction(ray_direction_){}
    };

    //3. sphere　
    struct Sphere{
        double radius;
        Vec position, emission, color; //position: sphere center
        int refl; //reflection type: 0=diffuse,1=specular,2=refractive

        //constructor
        Sphere(double radius_, Vec position_, Vec emission_, Vec color_, int refl): 
            radius(radius_),position(position_),emission(emission_),color(color_),refl(refl_){}

        //intersect(returns distance, 0=nohit)
        double intersect(const Ray &r) const{
            //solve t^2*d.d+2*t*(o-p).d+(o-p).(o-p)-R^2=0
            Vec op = position - r.ray_origin;
            double t, eps = 1e-4;
            double b = op.dot(r.ray_direction);
            double det = b*b-op.dot(op)+radius*radius;
            if(det<0) return 0;
            else det = sqrt(det);
            return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
        }
    }

    //4. scene:radius, position, emission, color, material(0=diffuse,1=specular,2=refractive)
    Sphere spheres[]={
        Sphere(1e5, Vec(1e5+1, 40.8, 81.6), Vec(),Vec(.75,.25,.25),0), //left wall
        Sphere(1e5, Vec(-1e5+99, 40.8, 81.6), Vec(),Vec(.25,.25,.75),0), //right wall
        Sphere(1e5, Vec(50, 40.8, 1e5), Vec(),Vec(.75,.75,.75),0), //back wall
        Sphere(1e5, Vec(50, 40.8, -1e5+170), Vec(),Vec(),0), //front wall
        Sphere(1e5, Vec(50, 1e5, 81.6), Vec(),Vec(.75,.75,.75),0), //bottom wall
        Sphere(1e5, Vec(50, -1e5+81.6, 81.6), Vec(),Vec(.75,.75,.75),0), //top wall
        Sphere(16.5, Vec(27, 16.5, 47), Vec(),Vec(1,1,1)*999,1), //mirror
        Sphere(16.5, Vec(73, 16.5, 78), Vec(),Vec(1,1,1)*999,2), //glass
        Sphere(1.5, Vec(50, 81.6-16.5, 81.6), Vec(4,4,4)*100,Vec(),0), //left wall
    };
    int numSpheres = sizeof(spheres)/sizeof(Sphere);

    //5. convert colors to displayable range [0,255]
    inline double clamp(double x){return x<0 ? 0 : x>1 ? 1 : x;}
    //convert float to int; gamma correction of 2.2
    inline int toInt(double x){ return int(pow(clamp(x), 1/2.2)*255+0.5);}

    //6. intersects ray with scene
    inline bool intersect(const Ray &r, double &t, int &id){
        double n = sizeof(spheres)/sizeof(Sphere);
        double d;
        double inf = t = 1e20;

        for(int i = int(n); i--;){
            if((d=spheres[i].intersect(r))&&d<t){
                t=d; //最近交点距离
                id=i; //球的id
            }
        }
        return t<inf; //true or false
    }

}
/*
namespace VCX::Labs::GettingStarted {
    static constexpr auto c_Sizes = std::to_array<std::pair<std::uint32_t, std::uint32_t>>({
        { 320U, 320U },
        { 640U, 640U } 
    });

    static constexpr auto c_SizeItems = std::array<char const *, 2> {
        "Small (320 x 320)",
        "Large (640 x 640)"
    };

    static constexpr auto c_BgItems = std::array<char const *, 3> {
        "White",
        "Black",
        "Checkboard"
    };

    CaseFixed::CaseFixed() :
        _textures(Engine::make_array<Engine::GL::UniqueTexture2D, c_Sizes.size()>(
            Engine::GL::SamplerOptions {
                .MinFilter = Engine::GL::FilterMode::Linear,
                .MagFilter = Engine::GL::FilterMode::Nearest
            })),
        _empty({
            Common::CreatePureImageRGB(c_Sizes[0].first, c_Sizes[0].second, { 2.f / 17, 2.f / 17, 2.f / 17 }),
            Common::CreatePureImageRGB(c_Sizes[1].first, c_Sizes[1].second, { 2.f / 17, 2.f / 17, 2.f / 17 })
        }) {
    }

    //控制边栏
    void CaseFixed::OnSetupPropsUI() {
        ImGui::Checkbox("Zoom Tooltip", &_enableZoom);
        //_recompute |= ImGui::Combo("Size", &_sizeId, c_SizeItems.data(), c_SizeItems.size());
        //_recompute |= ImGui::Combo("Background", &_bgId, c_BgItems.data(), c_BgItems.size());
    }

    //渲染
    Common::CaseRenderResult CaseFixed::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        auto const width = c_Sizes[_sizeId].first;
        auto const height = c_Sizes[_sizeId].second;
        if (_recompute) {
            _recompute = false;
            _task.Emplace([this, width, height]() {
                Common::ImageRGB image(0, 0);
                switch (_bgId) {
                case 0:
                    image = Common::CreatePureImageRGB(width, height, { 1., 1., 1. });
                    break;
                case 1:
                    image = Common::CreatePureImageRGB(width, height, { 0., 0., 0. });
                    break;
                case 2:
                    image = Common::CreateCheckboardImageRGB(width, height);
                    break;
                }
                for (std::size_t x = 0; x < width; ++x)
                    for (std::size_t y = 0; y < height; ++y)
                        if (x + y < width) image.At(x, y) = { 1., 0., 0. };
                return image;
            });
        }
        _textures[_sizeId].Update(_task.ValueOr(_empty[_sizeId]));
        return Common::CaseRenderResult {
            .Fixed     = true,
            .Image     = _textures[_sizeId],
            .ImageSize = c_Sizes[_sizeId],
        };
    }

    //交互
    void CaseFixed::OnProcessInput(ImVec2 const & pos) {
        auto         window  = ImGui::GetCurrentWindow();
        bool         hovered = false;
        bool         anyHeld = false;
        ImVec2 const delta   = ImGui::GetIO().MouseDelta;
        ImGui::ButtonBehavior(window->Rect(), window->GetID("##io"), &hovered, &anyHeld);
        if (! hovered) return;
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left) && delta.x != 0.f)
            ImGui::SetScrollX(window, window->Scroll.x - delta.x);
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left) && delta.y != 0.f)
            ImGui::SetScrollY(window, window->Scroll.y - delta.y);
        if (_enableZoom && ! anyHeld && ImGui::IsItemHovered())
            Common::ImGuiHelper::ZoomTooltip(_textures[_sizeId], c_Sizes[_sizeId], pos);
    }
}
*/