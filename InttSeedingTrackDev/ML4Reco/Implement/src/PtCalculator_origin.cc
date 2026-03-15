#include "PtCalculator.h"
#include <cmath>
#include <stdexcept>
#include <limits>
#include <array>
#include <memory>
#include <iostream>

#include <onnxruntime_cxx_api.h>   // ONNX Runtime C++ API
#include <fstream>
#include <nlohmann/json.hpp>       // Choose any json library (or the one you already have)

namespace {

// Convert 7 raw inputs -> 3 engineered features:
// raw order (as in PtCalculator.h comment):
//   [0] R34, [1] Z34, [2] R56, [3] Z56, [4] Rcalo, [5] Zcalo, [6] E
static std::vector<float> MLEprojRaw7_to_Feature3(const std::vector<float>& raw7)
{
    constexpr double eps = 1e-12;

    const double R0 = raw7.at(0), Z0 = raw7.at(1);
    const double R1 = raw7.at(2), Z1 = raw7.at(3);
    const double R2 = raw7.at(4), Z2 = raw7.at(5);
    const double E  = raw7.at(6);

    // direction eta from endpoint (0 -> 2) in RZ
    const double dR = (R2 - R0);
    const double dZ = (Z2 - Z0);
    const double vT = std::fabs(dR);              // RZ 里把 transverse 当作 |dR|
    const double eta_dir = std::asinh(dZ / (vT + eps));  // 对应 data.py: asinh(vz/(vT+eps))

    // kink angle between segments (0->1) and (1->2) in RZ
    const double v01x = (R1 - R0), v01z = (Z1 - Z0);
    const double v12x = (R2 - R1), v12z = (Z2 - Z1);

    const double n01 = std::sqrt(v01x*v01x + v01z*v01z);
    const double n12 = std::sqrt(v12x*v12x + v12z*v12z);
    const double denom = (n01 * n12) + eps;

    double cos_kink = (v01x*v12x + v01z*v12z) / denom;
    if (cos_kink > 1.0) cos_kink = 1.0;
    if (cos_kink < -1.0) cos_kink = -1.0;
    const double kink = std::acos(cos_kink);

    // logE
    const double logE = std::log(E + 1e-6);       // 对应 data.py: log(E+1e-6)

    return { static_cast<float>(eta_dir),
             static_cast<float>(logE),
             static_cast<float>(kink) };
}

} // anonymous namespace


namespace SiCaloPt {

PtCalculator::PtCalculator(const PtCalculatorConfig& cfg)
: m_cfg(cfg) {}

void PtCalculator::setConfig(const PtCalculatorConfig& cfg) {
    m_cfg = cfg;
}

bool PtCalculator::init(std::string* err)
{
    // EMD ML Model loading
    if (m_cfg.mlEMD_model_path) 
    {
        std::string err_string;
        auto fn = MakeOnnxInfer(*m_cfg.mlEMD_model_path, &err_string);
        if (!fn) 
        { 
            if (err) *err = "ONNX load mlEMD failed: " + err_string; 
            return false; 
        }
        setMLEMDInfer(std::move(fn));
        // scaler (optional)
        if (m_cfg.mlEMD_scaler_json) 
        {
            std::vector<float> mean, scale;
            if (!LoadScalerJson(*m_cfg.mlEMD_scaler_json, mean, scale, &err_string)) 
            {
                if (err) *err = "Load mlEMD scaler failed: " + err_string; 
                return false;
            }
            setMLEMDStandardizer(std::move(mean), std::move(scale));
        }
    }

    // Eproj ML model loading
    if (m_cfg.mlEproj_model_path) 
    {
        std::string err_string;
        auto fn = MakeOnnxInfer(*m_cfg.mlEproj_model_path, &err_string);
        if (!fn)
        { 
            if (err) *err = "ONNX load mlEproj failed: " + err_string; 
            return false; 
        }
        setMLEprojInfer(std::move(fn));
        if (m_cfg.mlEproj_scaler_json) 
        {
            std::vector<float> mean, scale;
            if (!LoadScalerJson(*m_cfg.mlEproj_scaler_json, mean, scale, &err_string)) 
            {
                if (err) *err = "Load mlEproj scaler failed: " + err_string; 
                return false;
            }
            setMLEprojStandardizer(std::move(mean), std::move(scale));
        }
    }

    // Combined Gate ML model loading
    if (m_cfg.mlCombined_model_path) 
    {
        std::string err_string;
        auto fn = MakeOnnxInfer(*m_cfg.mlCombined_model_path, &err_string);
        if (!fn) 
        { 
            if (err) *err = "ONNX load mlCombined failed: " + err_string; 
            return false; 
        }
        setMLCombinedInfer(std::move(fn));
        if (m_cfg.mlCombined_scaler_json) 
        {
            std::vector<float> mean, scale;
            if (!LoadScalerJson(*m_cfg.mlCombined_scaler_json, mean, scale, &err_string)) 
            {
                if (err) *err = "Load mlCombined scaler failed: " + err_string; 
                return false;
            }
            setMLCombinedStandardizer(std::move(mean), std::move(scale));
        }
    }

    return true;
}


PtResult PtCalculator::ComputePt(Method method, const AnyInput& input) const 
{
    try 
    {
        switch (method) 
        {
            case Method::MethodEMD:
                return ComputeEMD(std::get<InputEMD>(input));
            case Method::MethodEproj:
                return ComputeEproj(std::get<InputEproj>(input));
            case Method::MethodMLEMD:
                return ComputeMLEMD(std::get<InputMLEMD>(input));
            case Method::MethodMLEproj:
                return ComputeMLEproj(std::get<InputMLEproj>(input));
            case Method::MethodMLCombined:
                return ComputeMLCombined(std::get<InputMLCombined>(input));
            default:
                return PtResult{.pt_reco = NAN, .ok = false, .err = "Unknown method"};
        }
    } 
    catch (const std::bad_variant_access&) 
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "Input type does not match method"};
    }
}

PtResult PtCalculator::ComputeEMD(const InputEMD& in) const 
{
    if (in.EMD_Angle == 0.f) 
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "EMD_Angle is zero"};
    }

    int scenario = getScenario(in.EMD_Angle);
    if (scenario == -999)
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "Invalid scenario"};
    }

    if (consider_eta_dependence_on_EMDcompute)
    {
        switch (scenario)
        {
            case -10:
                // negative charge + projected cluster
                m_par_Ceta = 0.199854 + (0.00971966*x_eta*x_eta) + (-0.0177071*x_eta*x_eta*x_eta*x_eta);
                m_par_Power = -1.0;
                break;

            case -11:
                // negative charge + inner face
                m_par_Ceta = 0.198211 + (0.013064*x_eta*x_eta) + (-0.009812*x_eta*x_eta*x_eta*x_eta);
                m_par_Power = -1.0;
                break;

            case -12:
                // negative charge + detail
                m_par_Ceta = 0.197232 + (0.014244*x_eta*x_eta) + (-0.0188948*x_eta*x_eta*x_eta*x_eta);
                m_par_Power = -1.0;
                break;

            case 10:
                // positive charge + projected cluster
                m_par_Ceta = 0.203717 + (-0.000515*x_eta*x_eta) + (-0.00131*x_eta*x_eta*x_eta*x_eta);
                m_par_Power = -1.0;
                break;

            case 11:
                // positive charge + inner face
                m_par_Ceta = 0.197297 + (0.01297*x_eta*x_eta) + (-0.00979*x_eta*x_eta*x_eta*x_eta);
                m_par_Power = -1.0;
                break;

            case 12:
                // positive charge + detail
                m_par_Ceta = 0.200697 + (-0.0040371*x_eta*x_eta) + (-0.000426*x_eta*x_eta*x_eta*x_eta);
                m_par_Power = -1.0;
                break;

            default:
                return PtResult{.pt_reco = NAN, .ok = false, .err = "Unknown scenario"};
        }
    }
    const float pt = m_par_Ceta * std::pow(std::fabs(in.EMD_Angle), m_par_Power);
    return PtResult{.pt_reco = pt, .ok = true, .err = ""};
}


PtResult PtCalculator::ComputeEproj(const InputEproj& in) const 
{
    const float Distance_Z = std::fabs(in.Z_Calo - in.Z_vertex);
    const float Distance_R = std::fabs(in.Radius_Calo - in.Radius_vertex);
    const float Distance   = std::sqrt((Distance_R*Distance_R) + (Distance_Z*Distance_Z));

    if (Distance == 0.f) 
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "Distance is zero"};
    }
    const float pt = in.Energy_Calo*(Distance_R/Distance);
    return PtResult{.pt_reco = pt, .ok = true, .err = ""};
}

PtResult PtCalculator::ComputeMLEMD(const InputMLEMD& in) const 
{
    if (!m_mlEMD_infer) 
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEMD infer function not set"};
    }

    std::vector<float> x = in.features;

    if (x.size() != 2)
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEMD expects 2 features: {dphi, eta}"};
    }

    const float dphi = x[0];
    if (!std::isfinite(dphi) || dphi == 0.f)
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEMD dphi is invalid or zero"};
    }

    // convert external input dphi -> internal model input 1/dphi
    x[0] = 1.f / dphi;

    if (!m_mlEMD_mean.empty() && !m_mlEMD_scale.empty()) 
    {
        if (m_mlEMD_mean.size() != x.size() || m_mlEMD_scale.size() != x.size()) 
        {
            return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEMD standardizer dim mismatch"};
        }
        applyStandardize(x, m_mlEMD_mean, m_mlEMD_scale);
    }

    const float pt = m_mlEMD_infer(x);
    return PtResult{.pt_reco = pt, .ok = true, .err = ""};
}

PtResult PtCalculator::ComputeMLEproj(const InputMLEproj& in) const 
{
    if (!m_mlEproj_infer) 
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEproj infer function not set"};
    }

    std::vector<float> x = in.features;

    // --- compute pt0 (baseline) ---
    // training uses: pt0 = E * sin(theta),  pt_pred = pt0 * exp(delta)
    // If we only have R,Z, approximate direction in RZ plane:
    //   vT = |dR|, vmag = sqrt(dR^2 + dZ^2), sin(theta) ~ vT/vmag
    double pt0 = NAN;

    if (x.size() == 7)
    {
        const double R0 = x[0], Z0 = x[1];
        const double R2 = x[4], Z2 = x[5];
        const double E  = x[6];

        const double dR = (R2 - R0);
        const double dZ = (Z2 - Z0);

        const double vT   = std::fabs(dR);
        const double vmag = std::sqrt(dR*dR + dZ*dZ) + 1e-12;
        const double uT   = vT / vmag;   // ~ sin(theta)

        pt0 = E * uT;

        // 7 raw -> 3 engineered
        x = MLEprojRaw7_to_Feature3(x);
    }
    else if (x.size() == 3)
    {
        // x = [eta_dir, logE, kink]
        const double eta_dir = x[0];
        const double logE    = x[1];

        // E = exp(logE) - 1e-6  (inverse of log(E+1e-6))
        const double E = std::exp(logE) - 1e-6;

        // sin(theta) = 1/cosh(eta)  (since eta = asinh(vz/vT), equivalent to usual direction eta)
        const double uT = 1.0 / std::cosh(eta_dir);

        pt0 = E * uT;
    }
    else
    {
        return PtResult{.pt_reco = NAN, .ok = false,
                        .err = "MLEproj expects 7 raw inputs or 3 engineered features"};
    }

    if (!std::isfinite(pt0) || pt0 <= 0.0)
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEproj pt0 invalid"};
    }

    // --- standardize x (must match dim=3 after conversion) ---
    if (!m_mlEproj_mean.empty() && !m_mlEproj_scale.empty()) 
    {
        if (m_mlEproj_mean.size() != x.size() || m_mlEproj_scale.size() != x.size()) 
        {
            return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEproj standardizer dim mismatch"};
        }
        applyStandardize(x, m_mlEproj_mean, m_mlEproj_scale);
    }

    // --- infer delta, then pt = pt0 * exp(delta) ---
    const float delta_f = m_mlEproj_infer(x);
    const double delta  = static_cast<double>(delta_f);

    if (!std::isfinite(delta))
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEproj delta invalid"};
    }

    const double pt = pt0 * std::exp(delta);

    if (!std::isfinite(pt) || pt <= 0.0)
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLEproj pt invalid"};
    }

    return PtResult{.pt_reco = static_cast<float>(pt), .ok = true, .err = ""};
}

PtResult PtCalculator::ComputeMLCombined(const InputMLCombined& in) const 
{
    if (!m_mlCombined_infer) 
    {
        return PtResult{.pt_reco = NAN, .ok = false, .err = "MLCombined infer function not set"};
    }
    std::vector<float> x = in.features;
    if (!m_mlCombined_mean.empty() && !m_mlCombined_scale.empty()) 
    {
        if (m_mlCombined_mean.size() != x.size() || m_mlCombined_scale.size() != x.size()) 
        {
            return PtResult{.pt_reco = NAN, .ok = false, .err = "MLCombined standardizer dim mismatch"};
        }
        applyStandardize(x, m_mlCombined_mean, m_mlCombined_scale);
    }
    const float pt = m_mlCombined_infer(x);
    return PtResult{.pt_reco = pt, .ok = true, .err = ""};
}

void PtCalculator::setMLEMDInfer(InferFn fn) { m_mlEMD_infer = std::move(fn); }
void PtCalculator::setMLEprojInfer(InferFn fn) { m_mlEproj_infer = std::move(fn); }
void PtCalculator::setMLCombinedInfer(InferFn fn) { m_mlCombined_infer = std::move(fn); }

void PtCalculator::setMLEMDStandardizer(std::vector<float> mean, std::vector<float> scale) 
{
    m_mlEMD_mean  = std::move(mean);
    m_mlEMD_scale = std::move(scale);
}

void PtCalculator::setMLEprojStandardizer(std::vector<float> mean, std::vector<float> scale) 
{
    m_mlEproj_mean  = std::move(mean);
    m_mlEproj_scale = std::move(scale);
}

void PtCalculator::setMLCombinedStandardizer(std::vector<float> mean, std::vector<float> scale) 
{
    m_mlCombined_mean  = std::move(mean);
    m_mlCombined_scale = std::move(scale);
}

void PtCalculator::applyStandardize(std::vector<float>& x,
                                    const std::vector<float>& mean,
                                    const std::vector<float>& scale) 
{
    const size_t n = x.size();
    for (size_t i=0; i<n; ++i) 
    {
        if (scale[i] == 0.f) 
        {
            throw std::runtime_error("Scale value is zero during standardization");
        }
        x[i] = (x[i] - mean[i]) / scale[i];
    }
}

// --------------------------------------------------------------------
bool PtCalculator::LoadScalerJson(const std::string& path,
                                  std::vector<float>& mean,
                                  std::vector<float>& scale,
                                  std::string* err)
{
    try {
        std::ifstream fin(path);
        nlohmann::json js; fin >> js;
        if (!js.contains("mean") || !js.contains("scale")) {
            if (err) *err = "scaler json missing keys";
            return false;
        }
        mean  = js["mean"].get<std::vector<float>>();
        scale = js["scale"].get<std::vector<float>>();
        return true;
    } catch (const std::exception& e) {
        if (err) *err = e.what();
        return false;
    }
}

PtCalculator::InferFn PtCalculator::MakeOnnxInfer(const std::string& onnx_path,
                                                  std::string* err)
{
    try {
        static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ptcalc");
        Ort::SessionOptions so;
        so.SetIntraOpNumThreads(1);
        so.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

        auto sess = std::make_shared<Ort::Session>(env, onnx_path.c_str(), so);

        Ort::AllocatorWithDefaultOptions alloc;
        auto in_name  = std::string(sess->GetInputNameAllocated(0, alloc).get());
        auto out_name = std::string(sess->GetOutputNameAllocated(0, alloc).get());

        return [sess, in_name, out_name](const std::vector<float>& feats) -> float {
            const int64_t N = 1;
            const int64_t D = static_cast<int64_t>(feats.size());
            std::array<int64_t,2> shape{N, D};

            Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
            Ort::Value input = Ort::Value::CreateTensor<float>(
                mem,
                const_cast<float*>(feats.data()),
                feats.size(),
                shape.data(), shape.size()
            );

            const char* in_names[]  = { in_name.c_str()  };
            const char* out_names[] = { out_name.c_str() };

            auto outputs = sess->Run(Ort::RunOptions{nullptr}, in_names, &input, 1, out_names, 1);
            float* out_ptr = outputs.front().GetTensorMutableData<float>();
            return out_ptr[0];
        };
    } catch (const std::exception& e) {
        if (err) *err = e.what();
        return {};
    }
}

} // namespace SiCaloPt
