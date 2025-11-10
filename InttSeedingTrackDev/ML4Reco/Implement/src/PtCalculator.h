#ifndef PT_CALCULATOR_H
#define PT_CALCULATOR_H

#include <functional>
#include <string>
#include <variant>
#include <vector>
#include <optional>
#include <limits>
#include <array>
#include <memory>

namespace SiCaloPt {

// consistent return struct
struct PtResult 
{
    float pt_reco = std::numeric_limits<float>::quiet_NaN(); // reconstructed pt
    bool ok = false; // whether successful
    std::string err; // error message if any
};

// Five approaches enum
enum class Method 
{
    MethodEMD,   // EM Deflection Method
    MethodEproj, // Energy Projection Method
    MethodMLEMD, // ML EM Deflection Method
    MethodMLEproj, // ML Energy Projection Method
    MethodMLCombined // ML Combined Method
};

// Formula method input struct
struct InputEMD 
{
    float EMD_Angle = 0.f;
    float EMD_Eta = 0.f;
    float EMD_Radius = 0.f;
};

struct InputEproj 
{
    float Energy_Calo = 0.f;
    float Radius_Calo = 0.f;
    float Z_Calo = 0.f;
    float Radius_vertex = 0.f;
    float Z_vertex = 0.f;
};

// ML Model input struct, using vector is convenient for ONNX, length and order must be consistent with training
// Defination of length and order is in the tutorial file PtCalcMLTutorial.C
struct InputMLEMD 
{
    std::vector<float> features; 
};

struct InputMLEproj 
{
    std::vector<float> features;
};

struct InputMLCombined 
{
    std::vector<float> features;
};

// consistent input variant for all methods
using AnyInput = std::variant<InputEMD, InputEproj, InputMLEMD, InputMLEproj, InputMLCombined>;

// Config(optional) for ML: ML paths, standardizer params, etc.
struct PtCalculatorConfig 
{
    std::optional<std::string> mlEMD_model_path;
    std::optional<std::string> mlEproj_model_path;
    std::optional<std::string> mlCombined_model_path;

    std::optional<std::string> mlEMD_scaler_json;
    std::optional<std::string> mlEproj_scaler_json;
    std::optional<std::string> mlCombined_scaler_json;
};

// PtCalculator MAIN CLASS 
class PtCalculator {
public:
    PtCalculator() = default;
    explicit PtCalculator(const PtCalculatorConfig& cfg);

    // setting/configuration (does not auto-load external resources)
    void setConfig(const PtCalculatorConfig& cfg);

    // initialize (load models and scalers for ML methods)
    bool init(std::string* err = nullptr);

    // general func. for Pt calculation, dispatching according to method
    PtResult ComputePt(Method method, const AnyInput& input) const;

    // independent interfaces for finer control
    PtResult ComputeEMD(const InputEMD& in) const;
    PtResult ComputeEproj(const InputEproj& in) const;
    PtResult ComputeMLEMD(const InputMLEMD& in) const;
    PtResult ComputeMLEproj(const InputMLEproj& in) const;
    PtResult ComputeMLCombined(const InputMLCombined& in) const;

    // load ML infer, output is pt
    using InferFn = std::function<float(const std::vector<float>&)>;
    void setMLEMDInfer(InferFn fn);
    void setMLEprojInfer(InferFn fn);
    void setMLCombinedInfer(InferFn fn);

    //  optional: set standardizer for each ML model, or do standardization inside InferFn
    void setMLEMDStandardizer(std::vector<float> mean, std::vector<float> scale);
    void setMLEprojStandardizer(std::vector<float> mean, std::vector<float> scale);
    void setMLCombinedStandardizer(std::vector<float> mean, std::vector<float> scale);

    // EMD formula parameters setters
    void setParCeta(float v)  { m_par_Ceta = v; }
    void setParPower(float v) { m_par_Power = v; }

private:
    static void applyStandardize(std::vector<float>& x,
                                 const std::vector<float>& mean,
                                 const std::vector<float>& scale);

    static bool LoadScalerJson(const std::string& path,
                                std::vector<float>& mean,
                                std::vector<float>& scale,
                                std::string* err = nullptr);

    static SiCaloPt::PtCalculator::InferFn MakeOnnxInfer(const std::string& onnx_path, std::string* err = nullptr);

private:
    PtCalculatorConfig m_cfg;

    // ML Model infer functions
    InferFn m_mlEMD_infer;
    InferFn m_mlEproj_infer;
    InferFn m_mlCombined_infer;

    // optional: ML Scaler params for standardization, not exist for EMD and Combined model now
    std::vector<float> m_mlEMD_mean, m_mlEMD_scale;
    std::vector<float> m_mlEproj_mean, m_mlEproj_scale; 
    std::vector<float> m_mlCombined_mean, m_mlCombined_scale;

    // EMD formula parameters
    float m_par_Ceta = 0.2;
    float m_par_Power = -1.0;
};

} // namespace SiCaloPt

#endif // PT_CALCULATOR_H
