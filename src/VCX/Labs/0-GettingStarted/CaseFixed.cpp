#include <algorithm>
#include <array>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "Labs/0-GettingStarted/CaseFixed.h"
#include "Labs/Common/ImGuiHelper.h"

namespace VCX::Labs::GettingStarted {
    static constexpr auto c_Sizes = std::to_array<std::pair<std::uint32_t, std::uint32_t>>({
        { 512U, 384U },
        { 1024U, 768U } 
    });

    static constexpr auto c_SizeItems = std::array<char const *, 2> {
        "Small (512 x 384)",
        "Large (1024 x 768)"
    };

    static constexpr auto c_BgItems = std::array<char const *, 3> {
        "White",
        "Black",
        "Checkboard"
    };

    // 初始化可用采样数量的滑块范围
    static constexpr int c_SampsMin = 1;
    static constexpr int c_SampsMax = 100;

    CaseFixed::CaseFixed() :
        _textures(Engine::make_array<Engine::GL::UniqueTexture2D, c_Sizes.size()>(
            Engine::GL::SamplerOptions {
                .MinFilter = Engine::GL::FilterMode::Linear,
                .MagFilter = Engine::GL::FilterMode::Nearest
            })),
        _empty({
            Common::CreatePureImageRGB(c_Sizes[0].first, c_Sizes[0].second, { 2.f / 17, 2.f / 17, 2.f / 17 }),
            Common::CreatePureImageRGB(c_Sizes[1].first, c_Sizes[1].second, { 2.f / 17, 2.f / 17, 2.f / 17 })
        }) ,
         _samps(1), // 初始化采样值
        _recompute(true) // 默认需要重新计算
        {
    }

    //进度条：设置进度回调
    void CaseFixed::SetProgressCallback(ProgressCallback callback){
        _progressCallback = std::move(callback);
    }

    //控制边栏
    void CaseFixed::OnSetupPropsUI() {
        ImGui::Checkbox("Zoom Tooltip", &_enableZoom);
        _recompute |= ImGui::Combo("Size", &_sizeId, c_SizeItems.data(), c_SizeItems.size());
        //_recompute |= ImGui::Combo("Background", &_bgId, c_BgItems.data(), c_BgItems.size());
        
        // 添加滑块用于调整采样值
        int prevSamps = _samps;
        ImGui::SliderInt("Samples", &_samps, c_SampsMin, c_SampsMax);
        if (prevSamps != _samps) {
            _recompute = true; // 如果采样值发生变化，则标记需要重新计算
        }

        // 添加渲染进度条
        if (_isRendering) {
            ImGui::Text("Rendering...");
            ImGui::ProgressBar(_progress, ImVec2(-1, 0));
        } else {
            ImGui::Text("Ready");
    }
    }
 
    //渲染
    Common::CaseRenderResult CaseFixed::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        auto const width = 512;
        auto const height = 384;

        if(_recompute){
            _recompute = false;
            _task.Emplace([this, width, height](){
                Common::ImageRGB image (width, height);

                //进度条
                _isRendering = true;
                _progress = 0.0f;

                Vec* result = PathTracing(width, height, _samps, [this](float progress){
                    _progress = progress;
                }); 

                for (std::size_t y = 0; y < height; ++y){
                    for (std::size_t x = 0; x < width; ++x){
                        int i = y * width + x;
                        image.At(x,y)={correct(result[i].x), correct(result[i].y), correct(result[i].z)};
                    }
                }//遍历pixels结束
                delete[] result;

                //进度条
                _isRendering = false;

                return image;
            }
            );
        }
        
        //Common::ImageRGB image = Common::CreatePureImageRGB(width, height, { 0., 0., 0. });
        //Vec* c = PathTracing(512, 384, 1);        
    
        
        _textures[_sizeId].Update(_task.ValueOr(_empty[_sizeId]));
        return Common::CaseRenderResult {
            .Fixed     = true,
            .Image     = _textures[_sizeId],
            .ImageSize = {512,384},
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