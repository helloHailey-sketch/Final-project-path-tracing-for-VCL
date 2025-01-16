#include "Labs/0-GettingStarted/CaseResizable.h"
#include "Labs/Common/ImageRGB.h"
#include "Labs/Common/ImGuiHelper.h"

namespace VCX::Labs::GettingStarted {
    static constexpr auto c_Sizes = std::to_array<std::pair<std::uint32_t, std::uint32_t>>({
        { 512U, 384U },
        { 768U, 576U } 
    });

    static constexpr auto c_SizeItems = std::array<char const *, 2> {
        "Small (512 x 384)",
        "Large (1024 x 768)"
    };

    // 初始化可用采样数量的滑块范围
    static constexpr int c_SampsMin = 1;
    static constexpr int c_SampsMax = 100;

    CaseResizable::CaseResizable() :
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


    void CaseResizable::OnSetupPropsUI() {
        ImGui::Checkbox("Zoom Tooltip", &_enableZoom);
        //调整尺寸
        _recompute |= ImGui::Combo("Size", &_sizeId, c_SizeItems.data(), c_SizeItems.size());
        
        // 添加滑块用于调整采样值
        int prevSamps = _samps;
        ImGui::SliderInt("SPP/4", &_samps, c_SampsMin, c_SampsMax);
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

    Common::CaseRenderResult CaseResizable::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
       //auto const width = 512;
        //auto const height = 384;
        std::uint32_t width = c_Sizes[_sizeId].first;   
        std::uint32_t height = c_Sizes[_sizeId].second;

        if(_recompute){
            _recompute = false;
            _task.Emplace([this, width, height](){
                Common::ImageRGB image (width, height);

                //进度条
                _isRendering = true;
                _progress = 0.0f;

                Vec* result = PathTracing_painterly(width, height, _samps, [this](float progress){
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
            .ImageSize = {width,height},
        };
    }

    void CaseResizable::OnProcessInput(ImVec2 const& pos) {
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
    static constexpr auto c_PositionData = std::to_array<glm::vec2>({
        { -0.5, -0.5 },
        {  0  ,  0.5 },
        {  0.5, -0.5 },
    });

    static constexpr auto c_ColorData = std::to_array<glm::vec3>({
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    });

    CaseResizable::CaseResizable() :
         _program(
            Engine::GL::UniqueProgram({
                Engine::GL::SharedShader("assets/shaders/triangle.vert"),
                Engine::GL::SharedShader("assets/shaders/triangle.frag") })),
        _mesh(
            Engine::GL::VertexLayout()
                .Add<glm::vec2>("position", Engine::GL::DrawFrequency::Static, 0)
                .Add<glm::vec3>("color", Engine::GL::DrawFrequency::Static, 1)) {
        _mesh.UpdateVertexBuffer("position", Engine::make_span_bytes<glm::vec2>(c_PositionData));
        _mesh.UpdateVertexBuffer("color", Engine::make_span_bytes<glm::vec3>(c_ColorData));
    }

    void CaseResizable::OnSetupPropsUI() {
        ImGui::Checkbox("Zoom Tooltip", &_enableZoom);
    }

    Common::CaseRenderResult CaseResizable::OnRender(std::pair<std::uint32_t, std::uint32_t> const desiredSize) {
        _frame.Resize(desiredSize);
        gl_using(_frame);
        _mesh.Draw({ _program.Use() });
        return Common::CaseRenderResult {
            .Fixed     = false,
            .Flipped   = true,
            .Image     = _frame.GetColorAttachment(),
            .ImageSize = desiredSize,
        };
    }

    void CaseResizable::OnProcessInput(ImVec2 const& pos) {
        auto         window  = ImGui::GetCurrentWindow();
        bool         hovered = false;
        bool         anyHeld = false;
        ImGui::ButtonBehavior(window->Rect(), window->GetID("##io"), &hovered, &anyHeld);
        if (! hovered) return;
        if (_enableZoom && ! anyHeld && ImGui::IsItemHovered())
            Common::ImGuiHelper::ZoomTooltip(_frame.GetColorAttachment(), _frame.GetSize(), pos, true);
    }
}
*/
