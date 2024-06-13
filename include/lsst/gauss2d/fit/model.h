#ifndef LSST_GAUSS2D_FIT_MODEL_H
#define LSST_GAUSS2D_FIT_MODEL_H

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "lsst/gauss2d/evaluate.h"
#include "lsst/gauss2d/image.h"
#include "lsst/gauss2d/to_string.h"
#include "lsst/gauss2d/type_name.h"

#include "data.h"
#include "math.h"
#include "parameters.h"
#include "parametricmodel.h"
#include "param_filter.h"
#include "prior.h"
#include "psfmodel.h"
#include "source.h"
#include "util.h"

namespace lsst::gauss2d::fit {

static const std::string ERRMSG_PARAMS
        = "Were params fixed/freed since calling compute_*/evaluate?"
          "Or does your PSF model have free parameters?";

/// Valid modes to initialize lsst::gauss2d::GaussianEvaluator instances for.
enum class EvaluatorMode {
    image,          ///< Compute the model images only
    loglike,        ///< Compute the log likelihood only
    loglike_image,  ///< Compute the model images and log likelihood
    loglike_grad,   ///< Compute the gradients of the log likelihood
    jacobian        ///< Compute the model Jacobian
};

struct HessianOptions {
    /*
     * @param return_negative Whether the matrix should have all negative terms.
     *                        Should be set to true if the inverse Hessian is being used to estimate errors.
     * @param findiff_frac The value of the finite difference increment, as a fraction of the parameter value.
     * @param findiff_add The minimum of the finite difference increment (added to the fraction).
     */
    bool return_negative = true;
    double findiff_frac = 1e-4;
    double findiff_add = 1e-4;
};

typedef std::vector<std::shared_ptr<PsfModel>> PsfModels;
typedef std::vector<std::shared_ptr<Source>> Sources;
typedef std::vector<std::shared_ptr<Prior>> Priors;

/*
    A Model is a collection of Sources used to represent a model of the
    provided Data. The main purpose of Models is to evaluate themselves to
    generate output Images and/or gradients thereof to compare to the Data
    for fitting.

    Because Models are designed for repeated evaluations for fitting, one
    must explicitly set up the Evaluators to perform the exact type of
    evaluation required. In the future, a convenience "evaluate" function
    with arguments could be added, but it would in any case require
    storing and checking arguments against the previous call and therefore
    be less efficient.
*/
/**
 * @brief A collection of Sources comprising a ParametricModel of Data.
 *
 * @tparam T The type of the Image (usually float or double).
 * @tparam Image The class of image used in Data.
 * @tparam Indices The class of index map used in Data.
 * @tparam Mask The class of mask used in Data.
 */
template <typename T, typename Image, typename Indices, typename Mask>
class Model : public ParametricModel {
public:
    typedef GaussianEvaluator<T, Image, Indices> Evaluator;
    typedef std::map<std::reference_wrapper<const Channel>, std::vector<std::unique_ptr<const Gaussians>>>
            GaussiansMap;
    typedef Data<T, Image, Mask> ModelData;
    typedef typename ModelData::Observation Observation;

    /**
     * Construct a Model instance from Data, PsfModels, Sources and Priors.
     *
     * @param data The data to model.
     * @param psfmodels A vector of PSF models, ordered to match each Observation in data.
     * @param sources A vector of Source models.
     * @param priors A vector of Prior likelihoods.
     */
    explicit Model(std::shared_ptr<const ModelData> data, PsfModels& psfmodels, Sources& sources,
                   Priors& priors)
            : _data(std::move(data)),
              _size(_data == nullptr ? 0 : _data->size()),
              _size_priors(priors.size()) {
        if (_data == nullptr) throw std::invalid_argument("Model data can't be null");
        size_t size = this->size();
        if (psfmodels.size() != size) {
            throw std::invalid_argument("Model psfmodels.size()=" + std::to_string(psfmodels.size())
                                        + "!= data.size()=" + std::to_string(size));
        }

        _outputs.resize(size);
        _psfmodels.reserve(size);
        size_t i = 0;
        for (auto& psfmodel : psfmodels) {
            if (psfmodel == nullptr)
                throw std::invalid_argument("Model psfmodels[" + std::to_string(i) + "] can't be null");
            _psfcomps.push_back(psfmodel->get_gaussians());
            _psfmodels.push_back(psfmodel);
            i++;
        }

        _priors.reserve(_size_priors);
        i = 0;
        for (auto& prior : priors) {
            if (prior == nullptr)
                throw std::invalid_argument("Model priors[" + std::to_string(i) + "] can't be null");
            _priors.push_back(prior);
            i++;
        }

        _sources.reserve(sources.size());
        i = 0;
        for (auto& source : sources) {
            if (source == nullptr)
                throw std::invalid_argument("Model Sources[" + std::to_string(i) + "] can't be null");
            _sources.push_back(source);
            i++;
        }

        _gaussians_srcs = this->_get_gaussians_srcs();
        _factors_extra_in.resize(_size);
        _factors_grad_in.resize(_size);
        _map_extra_in.resize(_size);
        _map_grad_in.resize(_size);
    }

private:
    std::shared_ptr<const ModelData> _data;
    std::vector<std::unique_ptr<Evaluator>> _evaluators = {};
    ExtraParamFactors _factors_extra = {};
    std::vector<std::weak_ptr<Image>> _factors_extra_in;
    GradParamFactors _factors_grad = {};
    std::vector<std::weak_ptr<Image>> _factors_grad_in;
    GaussiansMap _gaussians_srcs;
    std::vector<std::shared_ptr<ImageArray<T, Image>>> _grads;
    bool _is_setup = false;
    std::vector<double> _likelihood_const_terms = {};
    std::vector<std::weak_ptr<Indices>> _map_extra_in;
    std::vector<std::weak_ptr<Indices>> _map_grad_in;
    EvaluatorMode _mode;
    size_t _n_params_free = 0;
    ParameterMap _offsets_params = {};
    std::vector<std::shared_ptr<Image>> _outputs = {};
    std::vector<std::shared_ptr<Prior>> _priors = {};
    std::vector<std::unique_ptr<const lsst::gauss2d::Gaussians>> _psfcomps;
    PsfModels _psfmodels;
    std::vector<std::shared_ptr<Image>> _outputs_prior = {};
    std::shared_ptr<Image> _residuals_prior = nullptr;
    size_t _size;
    size_t _size_priors;
    Sources _sources;

    void _check_obs_idx(size_t idx) const {
        if (!(idx < _size)) {
            throw std::invalid_argument("idx=" + std::to_string(idx)
                                        + " !< data->size()=" + std::to_string(_size));
        }
    }

    /// Return a map of source Gaussians for each Channel in Data
    GaussiansMap _get_gaussians_srcs() const {
        GaussiansMap gaussians_srcs{};
        size_t n_gaussians_src = 0;

        for (const Channel& channel : _data->get_channels()) {
            size_t n_gaussians_c = 0;
            gaussians_srcs[channel] = std::vector<std::unique_ptr<const Gaussians>>();
            auto& gaussians_src = gaussians_srcs[channel];
            for (const auto& src : this->get_sources()) {
                gaussians_src.push_back(src->get_gaussians(channel));
                n_gaussians_c += gaussians_src.back()->size();
            }
            if (n_gaussians_src == 0) {
                n_gaussians_src = n_gaussians_c;
            } else if (n_gaussians_src != n_gaussians_c) {
                throw std::logic_error("n_gaussians[" + channel.str() + "]=" + std::to_string(n_gaussians_c)
                                       + "!=n_gaussians=" + std::to_string(n_gaussians_src));
            }
        }
        return gaussians_srcs;
    }

    /**
     * Evaluate a single observation using the current EvaluatorMode in _mode
     *
     * @param idx The index of the observation.
     * @param print Whether to print diagnostic statements to stdout.
     * @return The log likelihood of the model for this observation
     */
    double _evaluate_observation(size_t idx, bool print = false, bool skip_zeroing = false) {
        if (print) std::cout << "Evaluating observation[" << idx << "]" << std::endl;
        if (_mode == EvaluatorMode::jacobian || _mode == EvaluatorMode::loglike_grad) {
            const auto& channel = _data->at(idx).get().get_channel();
            const size_t n_gaussians_psf = _psfcomps[idx]->size();
            const auto& gaussians_src = _gaussians_srcs.at(channel);

            const auto& factors_extra_in = _factors_extra_in[idx].lock();
            const auto& factors_grad_in = _factors_grad_in[idx].lock();

            size_t offset = 0;
            for (size_t i_psf = 0; i_psf < n_gaussians_psf; ++i_psf) {
                size_t i_src = 0;
                for (const auto& src : this->get_sources()) {
                    if (print)
                        std::cout << "Setting factors for psf[" << i_psf << "], src[" << i_src << "]"
                                  << std::endl;
                    const size_t n_gauss_src = gaussians_src[i_src++]->size();
                    if (n_gauss_src != src->get_n_gaussians(channel)) {
                        throw std::logic_error(
                                this->str() + "._data[" + channel.str() + "].get_n_gaussians(channel)="
                                + std::to_string(src->get_n_gaussians(channel))
                                + " != get_gaussians(channel).size()=" + std::to_string(n_gauss_src));
                    }
                    src->set_grad_param_factors(channel, _factors_grad, offset);
                    src->set_extra_param_factors(channel, _factors_extra, offset);
                    // Copy values to the actual input arrays used by the evaluator
                    const size_t row_max = offset + n_gauss_src;
                    for (size_t row = offset; row < row_max; ++row) {
                        if (print) {
                            std::cout << "factors_extra[" << row << "]=[";
                        }
                        for (size_t col = 0; col < lsst::gauss2d::N_EXTRA_FACTOR; ++col) {
                            factors_extra_in->set_value_unchecked(row, col, _factors_extra[row][col]);
                            if (print) {
                                std::cout << factors_extra_in->get_value_unchecked(row, col) << ",";
                            }
                        }
                        if (print) {
                            std::cout << "]\nfactors_grad[" << row << "]=[";
                        }
                        for (size_t col = 0; col < lsst::gauss2d::N_PARAMS_GAUSS2D; ++col) {
                            factors_grad_in->set_value_unchecked(row, col, _factors_grad[row][col]);
                            if (print) {
                                std::cout << factors_grad_in->get_value_unchecked(row, col) << ",";
                            }
                        }
                        if (print) {
                            std::cout << "]\n";
                        }
                    }
                    offset = row_max;
                }
            }
            if (!skip_zeroing) {
                auto& grads = *(this->_grads[idx]);
                for (auto& grad : grads) grad->fill(0);
            }
        }
        return _evaluators[idx]->loglike_pixel();
    }

    /**
     * Evaluate the prior using the current EvaluatorMode in _mode
     *
     * @param print Whether to print diagnostic statements to stdout.
     * @return The sum of the log likelihoods of the priors
     */
    double _evaluate_priors(bool print = false, bool normalize_loglike = false) {
        bool is_jacobian = _mode == EvaluatorMode::jacobian;
        double loglike = 0;
        size_t idx_resid = 0;

        if (print) std::cout << "Evaluating priors" << std::endl;
        for (const auto& prior : _priors) {
            auto eval = prior->evaluate(is_jacobian, normalize_loglike);
            loglike += eval.loglike;
            if (_residuals_prior != nullptr) {
                size_t idx_jac = idx_resid;
                size_t n_resid = eval.residuals.size();
                for (const double value : eval.residuals) {
                    this->_residuals_prior->set_value(0, idx_resid++, value);
                }
                // Reset the index back to the original start
                idx_resid = idx_jac;
                if (is_jacobian) {
                    for (const auto& [param, values] : eval.jacobians) {
                        if (values.size() != n_resid) {
                            throw std::logic_error("Prior=" + prior->str() + " param=" + param.get().str()
                                                   + " has n_values=" + std::to_string(values.size())
                                                   + " != residuals.size=" + std::to_string(n_resid));
                        }
                        const size_t idx_param = _offsets_params.at(param);
                        const auto& img = _outputs_prior.at(idx_param);
                        for (const double value : values) {
                            img->set_value(0, idx_jac++, value);
                        }
                        idx_jac = idx_resid;
                    }
                }
                idx_resid += n_resid;
            }
            if (print)
                std::cout << "Evaluated prior with loglike=" << eval.loglike << "; at idx_resid=" << idx_resid
                          << std::endl;
        }
        return loglike;
    }

    /**
     * Setup GaussianEvaluator instances for every Observation in Data using one EvaluatorMode.
     *
     * @param idx_obs The number index of the Observation in Data.
     * @param mode The EvaluatorMode for Evaluator instances.
     * @param outputs The list of Image instances to output to, ordered as Data.
     * @param residual The residual Image for the Evaluator to use.
     * @param print Whether to print diagnostic information to stdout.
     * @return A pair of the Evaluator instance, along with the associated output Image
     *         and gradient ImageArray.
     *
     * @note Different modes require different lengths of outputs. See setup_evaluators for details.
     */
    std::pair<std::unique_ptr<Evaluator>,
              std::pair<std::shared_ptr<Image>, std::shared_ptr<ImageArray<T, Image>>>>
    _make_evaluator(const size_t idx_obs, EvaluatorMode mode = EvaluatorMode::image,
                    std::vector<std::vector<std::shared_ptr<Image>>> outputs = {},
                    std::shared_ptr<Image> residual = nullptr, bool print = false) {
        _check_obs_idx(idx_obs);

        const Observation& observation = _data->at(idx_obs).get();
        const Channel& channel = observation.get_channel();
        const auto& psfcomps = _psfcomps.at(idx_obs);
        const size_t n_c = observation.get_n_cols();
        const size_t n_r = observation.get_n_rows();

        const size_t n_outputs = outputs.size();
        const bool has_outputs = n_outputs > 0;
        const auto* outputs_obs = has_outputs ? &(outputs.at(idx_obs)) : nullptr;

        const bool is_mode_image = mode == EvaluatorMode::image;
        const bool is_mode_loglike_grad = mode == EvaluatorMode::loglike_grad;
        const bool is_mode_jacobian = mode == EvaluatorMode::jacobian;
        const bool do_output = is_mode_image || (mode == EvaluatorMode::loglike_image);

        if (has_outputs) {
            if (is_mode_jacobian || is_mode_loglike_grad) {
                const size_t size_out = outputs_obs->size();
                size_t n_obsect = is_mode_jacobian ? (_n_params_free + 1) : this->size();
                if (size_out != n_obsect) {
                    throw std::invalid_argument("outputs[" + std::to_string(idx_obs)
                                                + "].size()=" + std::to_string(size_out) + "!=n_"
                                                + (is_mode_jacobian ? "free+1" : "obs") + "="
                                                + std::to_string(n_obsect));
                }
            }
            if (!(n_outputs > idx_obs)) {
                throw std::invalid_argument("outputs.size()=" + std::to_string(n_outputs)
                                            + " !> idx_obs=" + std::to_string(idx_obs));
            }
        }

        std::shared_ptr<const CoordinateSystem> coordsys = observation.get_image().get_coordsys_ptr_const();
        std::shared_ptr<Image> output
                = !do_output ? nullptr
                             : (has_outputs ? outputs_obs->at(0)
                                            : std::make_shared<Image>(n_r, n_c, nullptr, coordsys));

        auto data_eval = !is_mode_image ? observation.get_image_ptr_const() : nullptr;
        auto sigma_inv = !is_mode_image ? observation.get_sigma_inverse_ptr_const() : nullptr;

        const auto& gaussians_src = _gaussians_srcs.at(channel);
        const size_t n_gaussians_psf = psfcomps->size();

        for (size_t p = 0; p < n_gaussians_psf; ++p) {
            const auto& psfcomp = psfcomps->at(p);
            const auto value = psfcomp.get_integral_value();
            // TODO: Downgrade to warning and/or handle neatly
            if (!(value > 0)) {
                throw std::runtime_error("psfcomps[" + std::to_string(p) + "]=" + psfcomp.str()
                                         + " get_integral_value()=" + std::to_string(value) + "!>0");
            }
        }

        size_t n_gaussians_src = 0;
        for (const auto& gaussians : gaussians_src) n_gaussians_src += gaussians->size();
        const size_t n_gaussians_conv = n_gaussians_src * n_gaussians_psf;

        ConvolvedGaussians::Data data{};

        for (auto psfit = psfcomps->cbegin(); psfit < psfcomps->cend(); psfit++) {
            for (const auto& gaussians : gaussians_src) {
                for (const auto& gaussian : *gaussians) {
                    data.emplace_back(std::make_shared<ConvolvedGaussian>(gaussian, *psfit));
                }
            }
        }
        auto gaussians = std::make_shared<const ConvolvedGaussians>(data);
        if (gaussians->size() != n_gaussians_conv) {
            throw std::logic_error("gaussians.size()=" + std::to_string(gaussians->size())
                                   + "!= n_gaussians_conv=" + std::to_string(n_gaussians_conv)
                                   + "= n_gaussians_src=" + std::to_string(n_gaussians_src)
                                   + "= n_gaussians_psf=" + std::to_string(n_gaussians_psf)
                                   + "(is_mode_jacobian=" + std::to_string(is_mode_jacobian));
        }

        std::shared_ptr<const Indices> map_extra_in, map_grad_in;
        std::shared_ptr<const Image> factors_extra_in, factors_grad_in;

        if (is_mode_jacobian || is_mode_loglike_grad) {
            std::weak_ptr<Image> factors_extra_weak = _factors_extra_in[idx_obs];
            std::weak_ptr<Image> factors_grad_weak = _factors_grad_in[idx_obs];
            std::weak_ptr<Indices> map_extra_weak = _map_extra_in[idx_obs];
            std::weak_ptr<Indices> map_grad_weak = _map_grad_in[idx_obs];

            if (print) {
                this->_make_param_maps<true>(factors_extra_weak, factors_grad_weak, map_extra_weak,
                                             map_grad_weak, n_gaussians_conv, n_gaussians_psf, channel,
                                             coordsys, idx_obs, map_extra_in, map_grad_in, factors_extra_in,
                                             factors_grad_in);
            } else {
                this->_make_param_maps<false>(factors_extra_weak, factors_grad_weak, map_extra_weak,
                                              map_grad_weak, n_gaussians_conv, n_gaussians_psf, channel,
                                              coordsys, idx_obs, map_extra_in, map_grad_in, factors_extra_in,
                                              factors_grad_in);
            }
        }

        std::shared_ptr<ImageArray<T, Image>> grads;
        if (is_mode_jacobian) {
            typename ImageArray<T, Image>::Data arrays{};
            const size_t n_free = _n_params_free + 1;

            for (size_t idx_img = 0; idx_img < n_free; ++idx_img) {
                std::shared_ptr<Image> image;
                if (has_outputs) {
                    image = (*outputs_obs)[idx_img];
                    if (image == nullptr) {
                        throw std::invalid_argument("output[" + std::to_string(idx_img)
                                                    + "] can't be null for obs #" + std::to_string(idx_obs));
                    }
                    const auto& image_ref = observation.get_image();
                    std::string msg;
                    // We don't care about coordinate systems in this case
                    auto compatible = images_compatible(*image, image_ref, false, &msg);
                    if (!compatible) {
                        throw std::invalid_argument("outputs[" + std::to_string(idx_img) + "]=" + image->str()
                                                    + " incompatible with corresponding image="
                                                    + image_ref.str() + " for obs #" + std::to_string(idx_obs)
                                                    + " due to: " + msg);
                    }
                } else {
                    image = std::make_shared<Image>(n_r, n_c, nullptr, coordsys);
                }
                arrays.emplace_back(image);
            }
            grads = std::make_shared<ImageArray<T, Image>>(&arrays);
        } else if (is_mode_loglike_grad) {
            typename ImageArray<T, Image>::Data arrays{};
            if (print) {
                std::cout << "Setting grads[" << idx_obs
                          << "] to image of n_cols=_n_params_free + 1=" << _n_params_free + 1
                          << "; _offsets_params.size()=" << _offsets_params.size() << std::endl;
            }
            auto image = has_outputs ? outputs_obs->at(0)
                                     : std::make_shared<Image>(1, _n_params_free + 1, nullptr, coordsys);
            arrays.emplace_back(image);
            grads = std::make_shared<ImageArray<T, Image>>(&arrays);
        } else {
            grads = nullptr;
        }
        if (print) {
            std::cout << observation.str() << "\n";
            std::cout << gaussians->str() << "\n";
            std::cout << "output is null: " << (output == nullptr) << "\n";
        }

        // Ensure grads & output aren't nullptr after move
        auto grads2 = grads;
        auto output2 = output;
        std::pair<std::shared_ptr<Image>, std::shared_ptr<ImageArray<T, Image>>> second
                = {std::move(output2), std::move(grads2)};

        auto test = Evaluator(gaussians, data_eval, sigma_inv, output, residual, grads, map_grad_in,
                              factors_grad_in, map_extra_in, factors_extra_in);
        auto evaluator
                = std::make_unique<Evaluator>(gaussians, data_eval, sigma_inv, output, residual, grads,
                                              map_grad_in, factors_grad_in, map_extra_in, factors_extra_in
                                              //, background
                );
        return {std::move(evaluator), std::move(second)};
    }

    /// Make new factor/map inputs for GaussianEvaluator, if needed.
    template <bool print>
    void _make_param_maps(std::weak_ptr<Image> factors_extra_weak, std::weak_ptr<Image> factors_grad_weak,
                          std::weak_ptr<Indices> map_extra_weak, std::weak_ptr<Indices> map_grad_weak,
                          size_t n_gaussians_conv, size_t n_gaussians_psf, const Channel& channel,
                          std::shared_ptr<const CoordinateSystem> coordsys, size_t idx_obs,
                          std::shared_ptr<const Indices>& map_extra_in,
                          std::shared_ptr<const Indices>& map_grad_in,
                          std::shared_ptr<const Image>& factors_extra_in,
                          std::shared_ptr<const Image>& factors_grad_in) {
        unsigned int n_obsired = map_extra_weak.expired() + map_grad_weak.expired()
                                 + factors_grad_weak.expired() + factors_extra_weak.expired();
        const bool expired = n_obsired == 4;
        if (print) {
            std::cout << "expired=" << expired << "\n";
        }
        if (!((n_obsired == 0) || expired)) {
            throw std::logic_error("jacobian n_obsired=" + std::to_string(n_obsired) + " not in (0,4)");
        }
        auto map_extra_mut = expired ? std::make_shared<Indices>(n_gaussians_conv, 2, nullptr, coordsys)
                                     : map_extra_weak.lock();
        auto map_grad_mut = expired ? std::make_shared<Indices>(
                                    n_gaussians_conv, lsst::gauss2d::N_PARAMS_GAUSS2D, nullptr, coordsys)
                                    : map_grad_weak.lock();
        auto factors_extra_mut = expired ? std::make_shared<Image>(n_gaussians_conv, 3, nullptr, coordsys)
                                         : factors_extra_weak.lock();
        auto factors_grad_mut = expired ? std::make_shared<Image>(
                                        n_gaussians_conv, lsst::gauss2d::N_PARAMS_GAUSS2D, nullptr, coordsys)
                                        : factors_grad_weak.lock();
        ;

        ExtraParamMap map_extra = {};
        GradParamMap map_grad = {};
        ExtraParamFactors& factors_extra = _factors_extra;
        GradParamFactors& factors_grad = _factors_grad;
        factors_extra.clear();
        factors_grad.clear();

        map_extra.reserve(n_gaussians_conv);
        map_grad.reserve(n_gaussians_conv);
        factors_extra.reserve(n_gaussians_conv);
        factors_grad.reserve(n_gaussians_conv);

        /*
            Make (and set) the param map; only make the factors (re-evaluated each time)
            The maps depend on the order of free parameters.
            the factors depend on their values, and can change every evaluation.
        */
        for (size_t i_psf = 0; i_psf < n_gaussians_psf; ++i_psf) {
            for (const auto& src : this->get_sources()) {
                src->add_grad_param_map(channel, map_grad, _offsets_params);
                src->add_extra_param_map(channel, map_extra, map_grad, _offsets_params);

                src->add_grad_param_factors(channel, factors_grad);
                src->add_extra_param_factors(channel, factors_extra);
            }
            if constexpr (print) {
                std::cout << "_offsets_params set to:" << std::endl;
                std::cout << str_map_ref<true>(_offsets_params) << std::endl;
            }
        }
        if (map_grad.size() != n_gaussians_conv) {
            throw std::logic_error("map_grad.size()=" + std::to_string(map_grad.size())
                                   + "!= n_gaussians_conv=" + std::to_string(n_gaussians_conv));
        }
        /*
            Validate the map if there is already another Jacobian evaluator.

            Do not support multiple evaluators set up with different maps (yet).
        */
        if (!expired) {
            std::string sizes_mut = std::to_string(map_extra_mut->size()) + ","
                                    + std::to_string(map_grad_mut->size()) + ","
                                    + std::to_string(factors_extra_mut->size()) + ","
                                    + std::to_string(factors_grad_mut->size());

            std::string sizes = std::to_string(map_extra.size()) + "," + std::to_string(map_grad.size()) + ","
                                + std::to_string(factors_extra.size()) + ","
                                + std::to_string(factors_grad.size());

            if (sizes != sizes_mut) {
                throw std::runtime_error(
                    "make_evaluator Jacobian map_extra,map_grad,factors_extra,factors_grad sizes_new="
                    + sizes + "!=sizes_mut=" + sizes_mut + "; did you make a new evaluator with different"
                                                           "free parameters from an old one?");
            }

            std::array<std::string, 4> errmsgs = {"", "", "", ""};

            size_t idx_err = 0;
            if constexpr (print) std::cout << "map_extra=[" << std::endl;
            for (size_t idx_row = 0; idx_row < map_extra.size(); idx_row++) {
                const auto& row = map_extra[idx_row];
                for (size_t col = 0; col < lsst::gauss2d::N_EXTRA_MAP; col++) {
                    const auto value_mut = map_extra_mut->get_value_unchecked(idx_row, col);
                    if ((idx_err < 2) && (row[col] != value_mut))
                        errmsgs[idx_err++] = "map_extra[" + std::to_string(idx_row) + "]["
                                             + std::to_string(col) + "]=" + std::to_string(row[col])
                                             + " != old=" + std::to_string(value_mut) + "\n";
                }
                if constexpr (print) {
                    std::cout << "    ";
                    stream_iter_ref(row, std::cout);
                    std::cout << std::endl;
                }
            }
            if constexpr (print) {
                std::cout << "]" << std::endl << "map_grad=[" << std::endl;
            }
            for (size_t idx_row = 0; idx_row < map_grad.size(); idx_row++) {
                const auto& row = map_grad[idx_row];
                for (size_t col = 0; col < lsst::gauss2d::N_EXTRA_MAP; col++) {
                    const auto value_mut = map_grad_mut->get_value_unchecked(idx_row, col);
                    if ((idx_err < 4) && (row[col] != value_mut))
                        errmsgs[idx_err++] = "map_grad[" + std::to_string(idx_row) + "]["
                                             + std::to_string(col) + "]=" + std::to_string(row[col])
                                             + " != old=" + std::to_string(value_mut) + "\n";
                }
                if constexpr (print) {
                    std::cout << "    ";
                    stream_iter_ref(row, std::cout);
                    std::cout << std::endl;
                }
            }
            if constexpr (print) std::cout << "]" << std::endl;
            std::string errmsg_extra = errmsgs[0] + errmsgs[1];
            std::string errmsg_grad = errmsgs[0] + errmsgs[1];
            if (!errmsg_extra.empty() || !errmsg_grad.empty()) {
                std::string errmsg
                        = "make_evaluator Jacobian map_extra/map_grad have different values"
                          " from old; did you make a new evaluator with different free parameters from an "
                          " old one? Errors:\n"
                          + errmsg_extra;
                if (!errmsg_extra.empty() && !errmsg_grad.empty()) errmsg += "...\n";
                errmsg += errmsg_grad;
                throw std::runtime_error(errmsg);
            }
        } else {
            /*
                Maps are stored as vectors in gauss2dfit (to avoid unnecessary templating)
                The inputs to gauss2d's Evaluator are templated Images, so copy the values
            */
            if constexpr (print) {
                std::cout << "factors/maps: {\n";
            }
            for (size_t row = 0; row < n_gaussians_conv; ++row) {
                for (size_t col = 0; col < lsst::gauss2d::N_EXTRA_MAP; ++col) {
                    map_extra_mut->set_value_unchecked(row, col, (map_extra)[row][col]);
                }
                for (size_t col = 0; col < lsst::gauss2d::N_EXTRA_FACTOR; ++col) {
                    factors_extra_mut->set_value_unchecked(row, col, (factors_extra)[row][col]);
                }
                for (size_t col = 0; col < lsst::gauss2d::N_PARAMS_GAUSS2D; ++col) {
                    map_grad_mut->set_value_unchecked(row, col, (map_grad)[row][col]);
                    factors_grad_mut->set_value_unchecked(row, col, (factors_grad)[row][col]);
                }
                if constexpr (print) {
                    std::cout << "    row[" << row << "]: map_extra=[";
                    stream_iter_ref(map_extra[row], std::cout);
                    std::cout << ", factors_extra=";
                    stream_iter_ref(factors_extra[row], std::cout);
                    std::cout << ", map_grad=";
                    stream_iter_ref(map_grad[row], std::cout);
                    std::cout << ", factors_grad=";
                    stream_iter_ref(factors_grad[row], std::cout);
                    std::cout << "\n";
                }
            }
            if constexpr (print) std::cout << "}" << std::endl;
        }
        _factors_extra_in[idx_obs] = factors_extra_mut;
        _factors_grad_in[idx_obs] = factors_grad_mut;

        map_extra_in = std::move(map_extra_mut);
        map_grad_in = std::move(map_grad_mut);
        factors_extra_in = std::move(factors_extra_mut);
        factors_grad_in = std::move(factors_grad_mut);
    }

public:
    void add_extra_param_map(const Channel& channel, ExtraParamMap& map_extra, const GradParamMap& map_grad,
                             ParameterMap& offsets) const override {
        for (auto& source : _sources) source->add_extra_param_map(channel, map_extra, map_grad, offsets);
    }
    void add_extra_param_factors(const Channel& channel, ExtraParamFactors& factors) const override {
        for (auto& source : _sources) source->add_extra_param_factors(channel, factors);
    }
    void add_grad_param_map(const Channel& channel, GradParamMap& map, ParameterMap& offsets) const override {
        for (auto& source : _sources) source->add_grad_param_map(channel, map, offsets);
    }
    void add_grad_param_factors(const Channel& channel, GradParamFactors& factors) const override {
        for (auto& source : _sources) source->add_grad_param_factors(channel, factors);
    }

    /**
     * Compute the gradient (partial first derivative) of the log-likelihood for each free parameter.
     *
     * @param include_prior Whether to include the prior likelihood(s) in the gradients.
     * @param print Whether to print diagnostic/debugging information.
     * @param verify Whether to verify the values by comparing to finite differences.
     * @param findiff_frac The value of the finite difference increment, as a fraction of the parameter value.
     * @param findiff_add The minimum value of the finite difference increment.
     * @param rtol The allowed relative tolerance in the Jacobian as compared to the finite difference.
     * @param atol The allowed absolute tolerance in the Jacobian as compared to the finite difference.
     * @return The gradient of the log likelihoods for the free parameters of the model, in the order
     *         returned by get_parameters.
     */
    std::vector<double> compute_loglike_grad(bool include_prior = true, bool print = false,
                                             bool verify = false, double findiff_frac = 1e-5,
                                             double findiff_add = 1e-5, double rtol = 1e-3,
                                             double atol = 1e-8) {
        this->setup_evaluators(EvaluatorMode::loglike_grad, {}, {}, {}, nullptr, true, print);
        this->evaluate(print);

        auto filter_free = g2f::ParamFilter{false, true, true, true};
        auto params_free = this->get_parameters_const_new(&filter_free);
        params_free = nonconsecutive_unique<ParamBaseCRef>(params_free);

        const size_t n_params = this->_offsets_params.size();
        if (n_params != params_free.size()) {
            throw std::runtime_error("_offsets_params=" + str_map_ref<true>(_offsets_params)
                                     + ".size() = " + std::to_string(n_params)
                                     + " != params_free=" + str_iter_ref<true>(params_free) + ".size() + 1 ="
                                     + std::to_string(params_free.size()) + "; " + ERRMSG_PARAMS);
        }
        std::map<ParamBaseCRef, size_t> param_idx = {};
        for (size_t idx_param = 0; idx_param < n_params; idx_param++) {
            const auto param = params_free[idx_param];
            if (_offsets_params.find(param) == _offsets_params.end()) {
                throw std::runtime_error(param.get().str() + " not found in _offsets_params\n"
                                         + ERRMSG_PARAMS);
            }
            param_idx[param] = idx_param;
        }
        if (print) {
            std::cout << "params_free=" << str_iter_ref<true>(params_free) << std::endl;
            std::cout << "_offsets_params=" << str_map_ref<true>(_offsets_params) << std::endl;
            std::cout << "param_idx=" << str_map_ref<true>(param_idx) << std::endl;
        }
        std::vector<double> loglike_grads(n_params, 0.);

        size_t idx_obs = 0;
        for (const auto& grads_obs : _grads) {
            const auto& values = grads_obs->at(0);
            const size_t n_cols = values.get_n_cols();
            if (n_cols != n_params + 1) {
                throw std::logic_error(this->_data->at(idx_obs).get().str()
                                       + " loglike_grads n_cols=" + std::to_string(n_cols)
                                       + " != n_params (from offsets) = " + std::to_string(n_params + 1));
            }
            for (const auto& [paramref, col] : _offsets_params) {
                loglike_grads.at(param_idx[paramref]) += values.get_value(0, col);
            }
            idx_obs++;
        }
        if (include_prior) {
            for (const auto& prior : _priors) {
                auto result = prior->evaluate(true);

                // TODO: Confirm that this works with ShapePrior
                for (const auto& paramref : params_free) {
                    double dll_dx = result.compute_dloglike_dx(paramref);
                    loglike_grads.at(param_idx[paramref]) += dll_dx;
                }
            }
        }
        if (verify) {
            std::string errmsg;
            this->setup_evaluators(g2f::EvaluatorMode::loglike);
            auto loglike = this->evaluate();

            auto params = this->get_parameters_new(&filter_free);
            params = nonconsecutive_unique<ParamBaseRef>(params);

            for (size_t idx_param = 0; idx_param < n_params; idx_param++) {
                auto& param = params[idx_param].get();
                double value = param.get_value_transformed();
                double diff = std::copysign(std::abs(value * findiff_frac) + findiff_add,
                                            loglike_grads[idx_param]);
                diff = finite_difference_param(param, diff);
                std::vector<double> loglike_new_plus;
                std::vector<double> loglike_new_minus;
                try {
                    param.set_value_transformed(value + diff / 2.);
                    loglike_new_plus = this->evaluate();
                    param.set_value_transformed(value - diff / 2.);
                    loglike_new_minus = this->evaluate();
                } catch (std::runtime_error& e) {
                    param.set_value_transformed(value + diff);
                    loglike_new_plus = this->evaluate();
                    loglike_new_minus = loglike;
                }
                auto dloglike = (sum_iter(loglike_new_plus) - sum_iter(loglike_new_minus)) / diff;
                auto close = isclose(dloglike, loglike_grads[idx_param], rtol, atol);
                param.set_value_transformed(value);
                if (!close.isclose) {
                    double dll_dx_exact = loglike_grads[idx_param];
                    double dlp_dx_sum = 0;
                    for (const auto& prior : _priors) {
                        double dlp_dx = prior->evaluate(true).compute_dloglike_dx(param, true);
                        dlp_dx_sum += dlp_dx;
                    }
                    double dlp_dx_findiff = (loglike_new_plus.back() - loglike_new_minus.back()) / diff;
                    errmsg += param.str() + " failed loglike_grad verification; isclose=" + close.str()
                              + " from findiff=" + to_string_float(dloglike) + " vs "
                              + to_string_float(dll_dx_exact)
                              + " (diff = " + to_string_float(dll_dx_exact - dloglike)
                              + " , ratio = " + to_string_float(dll_dx_exact / dloglike)
                              + "); dll_dx_prior=" + to_string_float(dlp_dx_sum)
                              + " vs findiff: " + to_string_float(dlp_dx_findiff) + "\n";
                }
            }
            if (errmsg.size() > 0) throw std::logic_error(errmsg);
        }

        return loglike_grads;
    }

    /**
     * Compute the Hessian matrix (second order partial derivatives) of the log likehood.
     *
     * @param transformed Whether the matrix should be computed for transformed parameters or not
     *                    If not, parameter transforms are temporarily removed.
     * @param include_prior Whether to include the prior likelihood(s) in the Hessian.
     * @param options Options for computing the Hessian via finite differencing of loglikelihood gradients.
     *                If null, the Hessian is estimated as J^T J (where J is the Jacobian).
     * @return The Hessian matrix for all free parameters.
     */
    std::unique_ptr<Image> compute_hessian(bool transformed = false, bool include_prior = true,
                                           std::optional<HessianOptions> options = std::nullopt,
                                           bool print = false) {
        this->setup_evaluators(EvaluatorMode::loglike_grad);

        auto filter_free = g2f::ParamFilter{false, true, true, true};
        auto params_free = this->get_parameters_new(&filter_free);
        params_free = nonconsecutive_unique<ParamBaseRef>(params_free);
        const size_t n_params_free = params_free.size();

        std::vector<std::shared_ptr<const parameters::Transform<double>>> transforms = {};
        transforms.reserve(n_params_free);
        if (!transformed) {
            for (auto& paramref : params_free) {
                auto& param = paramref.get();
                transforms.emplace_back(param.get_transform_ptr());
                double value = param.get_value();
                param.set_transform(nullptr);
                if (param.get_value() != value) {
                    throw std::logic_error("Param " + param.str() + " changed value from "
                                           + std::to_string(value) + " after dropping transform");
                }
            }
        }

        auto hessian = std::make_unique<Image>(n_params_free, n_params_free);
        const auto n_obs = this->size();
        std::vector<size_t> param_idx(n_params_free);

        size_t n_integral = 0;
        size_t n_frac = 0;

        for (size_t idx_param = 0; idx_param < n_params_free; idx_param++) {
            const auto param = params_free[idx_param];
            const auto found = _offsets_params.find(param);
            if (found == _offsets_params.end()) {
                throw std::runtime_error(param.get().str() + " not found in _offsets_params\n"
                                         + ERRMSG_PARAMS);
            }
            param_idx[idx_param] = (*found).second;
            if (not options) {
                const auto name = param.get().get_name();
                n_integral += name == IntegralParameterD::_name;
                n_frac += name == ProperFractionParameterD::_name;
            }
        }

        if (options) {
            const HessianOptions& opt = *options;

            auto loglike_grad = this->compute_loglike_grad(include_prior, false, false);

            size_t idx_param = 0;
            for (auto& paramref : params_free) {
                auto& param = paramref.get();
                const double value = param.get_value_transformed();
                // Make the finite differencing go away from the possible peak likelihood
                // This is probably a good idea even if return_negative isn't necessary
                double diff = std::copysign(std::abs(value * opt.findiff_frac) + opt.findiff_add,
                                            -loglike_grad[idx_param]);
                diff = finite_difference_param(param, diff);
                double diffabs = std::abs(diff);
                this->evaluate();

                auto loglike_grad2 = this->compute_loglike_grad(include_prior, false, false);
                for (size_t idx_col = 0; idx_col < n_params_free; ++idx_col) {
                    double dll_dx = loglike_grad2[idx_col] - loglike_grad[idx_col];
                    dll_dx *= opt.return_negative ? (1 - 2 * (dll_dx > 0)) / diffabs : 1 / diff;
                    hessian->set_value_unchecked(idx_param, idx_col, dll_dx);
                }
                param.set_value_transformed(value);
                idx_param++;
            }
        } else {
            if ((n_integral > 0) && (n_frac > 0)) {
                throw std::runtime_error(
                        "Cannot compute_hessian without options for model with free Integral and "
                        "ProperFraction "
                        "parameters (this is a bug); please supply options to use finite differencing.");
            }

            this->setup_evaluators(EvaluatorMode::jacobian);
            this->evaluate();

            hessian->fill(0);

            for (size_t idx_obs = 0; idx_obs < n_obs; ++idx_obs) {
                const auto& jacs = this->_grads.at(idx_obs);
                for (size_t idx_param1 = 0; idx_param1 < n_params_free; ++idx_param1) {
                    const auto& jac1 = jacs->at(param_idx[idx_param1]);
                    const size_t n_rows = jac1.get_n_rows();
                    const size_t n_cols = jac1.get_n_cols();

                    for (size_t idx_param2 = 0; idx_param2 < n_params_free; ++idx_param2) {
                        const auto& jac2 = jacs->at(param_idx[idx_param2]);
                        const size_t n_rows2 = jac2.get_n_rows();
                        const size_t n_cols2 = jac2.get_n_cols();

                        if ((n_rows != n_rows2) || (n_cols != n_cols2)) {
                            throw std::logic_error("n_rows,n_cols=" + to_string_float(n_rows) + ","
                                                   + to_string_float(n_cols)
                                                   + "!= n_rows2,n_cols2=" + to_string_float(n_rows2) + ","
                                                   + to_string_float(n_cols2));
                        }
                        double dx = 0;
                        for (size_t row = 0; row < n_rows; ++row) {
                            for (size_t col = 0; col < n_cols; ++col) {
                                // Note the sign convention of negative diagonal values
                                // at maximum likelihood
                                dx -= jac1.get_value(row, col) * jac2.get_value(row, col);
                            }
                        }
                        if (print) {
                            std::cout << "calling hessian->add_value(" << idx_param1 << "," << idx_param2
                                      << "," << dx << ")" << std::endl;
                        }
                        hessian->add_value(idx_param1, idx_param2, dx);
                    }
                }
            }
            if (include_prior && (_priors.size() > 0)) {
                std::map<ParamBaseCRef, size_t> param_idx_map = {};
                for (size_t idx_param = 0; idx_param < n_params_free; idx_param++) {
                    const auto param = params_free[idx_param];
                    if (_offsets_params.find(param) == _offsets_params.end()) {
                        throw std::runtime_error(param.get().str() + " not found in _offsets_params\n"
                                                 + ERRMSG_PARAMS);
                    }
                    param_idx_map[param] = idx_param;
                }

                size_t prior_size = 0;
                for (const auto& prior : _priors) {
                    prior_size += prior->size();
                }
                auto dll_prior = std::make_unique<Image>(n_params_free, prior_size);
                dll_prior->fill(0);

                size_t idx_prior = 0;
                for (const auto& prior : _priors) {
                    auto result = prior->evaluate(true);

                    for (const auto& [param_cref, values_base] : result.jacobians) {
                        size_t idx_residual = 0;
                        size_t idx_param = param_idx_map[param_cref];
                        for (const auto& value_base : values_base) {
                            dll_prior->add_value(idx_param, idx_prior + idx_residual++, value_base);
                        }
                    }
                    idx_prior += prior->size();
                }

                for (size_t row = 0; row < n_params_free; ++row) {
                    for (size_t row2 = 0; row2 < n_params_free; ++row2) {
                        double dll_dx = 0;
                        for (size_t col = 0; col < prior_size; ++col) {
                            dll_dx -= dll_prior->get_value(row, col) * dll_prior->get_value(row2, col);
                        }
                        hessian->add_value(row, row2, dll_dx);
                        if (print) {
                            std::cout << "calling hessian->add_value(" << row << "," << row2 << "," << dll_dx
                                      << ") for prior" << std::endl;
                        }
                    }
                }
            }
        }

        if (!transformed) {
            size_t idx_param = 0;
            for (auto& paramref : params_free) {
                auto& param = paramref.get();
                double value = param.get_value();
                param.set_transform(transforms[idx_param]);
                if (param.get_value() != value) {
                    throw std::logic_error("Param " + param.str() + " changed value from "
                                           + std::to_string(value) + " after dropping transform");
                }
                idx_param++;
            }
        }

        return hessian;
    }

    /**
     * Evaluate the model for every Observation in _data.
     *
     * @param print Whether to print diagnostic statements to stdout.
     * @param normalize_loglike Whether to include the normalizing (variance-dependent) term in the log
     *                          likelihood. If false, the log likelihood for a model with no residuals is 0.
     * @return The log likelihood of each observation, followed by the summed log likelihood of the priors.
     */
    std::vector<double> evaluate(bool print = false, bool normalize_loglike = false) {
        if (!_is_setup) throw std::runtime_error("Can't call evaluate before setup_evaluators");

        std::vector<double> result(_size + 1);
        bool is_loglike_grad = this->_mode == EvaluatorMode::loglike_grad;

        if (is_loglike_grad) {
            // loglike_grad adds values, and they have to be reset to zero first
            // This must be done for all observations in this case, because
            // most parameters change grads in all observations
            for (size_t idx = 0; idx < _size; ++idx) {
                this->_grads[idx]->at(0).fill(0);
            }
        }

        for (size_t idx = 0; idx < _size; ++idx) {
            result[idx] = this->_evaluate_observation(idx, print, is_loglike_grad);
        }
        if ((this->_mode == EvaluatorMode::loglike) || is_loglike_grad
            || (this->_mode == EvaluatorMode::loglike_image) || (this->_mode == EvaluatorMode::jacobian)) {
            result[_size] = this->_evaluate_priors(print, normalize_loglike);
        }

        return result;
    }

    /**
     * Evaluate a single observation with the given index in _data.
     *
     * @param idx The numeric index of the observation.
     * @return The log likelihood of the model for that observation.
     */
    double evaluate_observation(size_t idx) {
        _check_obs_idx(idx);
        return _evaluate_observation(_evaluators[idx]);
    }

    /// Return _data
    std::shared_ptr<const ModelData> get_data() const { return _data; }

    EvaluatorMode get_mode() const { return _mode; }

    std::unique_ptr<const lsst::gauss2d::Gaussians> get_gaussians(const Channel& channel) const override {
        std::vector<std::optional<const lsst::gauss2d::Gaussians::Data>> in;
        // TODO: Rethink this; it's not sufficient unless sources are single component gaussians
        in.reserve(_sources.size());

        for (auto& source : _sources) {
            auto data = source->get_gaussians(channel)->get_data();
            in.emplace_back(data);
        }

        return std::make_unique<lsst::gauss2d::Gaussians>(in);
    }

    /**
     * Get the constant (variance-dependent) terms of the log likelihood for each observation.
     *
     * @return A vector with the constant terms for each observation, and a final term for the priors.
     */
    std::vector<double> get_loglike_const_terms() {
        const size_t n_data = this->size();
        if (this->_likelihood_const_terms.empty()) {
            this->_likelihood_const_terms.resize(n_data + 1);

            for (size_t idx = 0; idx < n_data; ++idx) {
                const auto& sigma_inv = this->_data[idx].get_sigma_inverse();
                double loglike = 0;
                const size_t n_rows = sigma_inv.get_n_rows();
                const size_t n_cols = sigma_inv.get_n_cols();
                for (size_t col = 0; col < n_cols; ++col) {
                    for (size_t row = 0; row < n_rows; ++row) {
                        loglike += LOG_1 - log(sigma_inv._get_value(row, col));
                    }
                }
                _likelihood_const_terms[idx] = loglike;
            }
        }
        double loglike = 0;
        for (const auto& prior : this->_priors) {
            auto values = prior->get_loglike_const_terms();
            for (const auto value : values) loglike += value;
        }
        _likelihood_const_terms[n_data] = loglike;
        return _likelihood_const_terms;
    }

    size_t get_n_gaussians(const Channel& channel) const override {
        size_t n_g = 0;
        for (auto& source : _sources) n_g += source->get_n_gaussians(channel);
        return n_g;
    }

    /// Return _outputs (output Image instances for each Observation in _data)
    std::vector<std::shared_ptr<Image>> get_outputs() const { return _outputs; }

    std::vector<std::pair<ParamBaseCRef, size_t>> get_offsets_parameters() const {
        std::vector<std::pair<ParamBaseCRef, size_t>> offsets = {};
        for (auto& [param, offset] : _offsets_params) {
            offsets.emplace_back(param, offset);
        }
        return offsets;
    }

    ParamRefs& get_parameters(ParamRefs& params, ParamFilter* filter = nullptr) const override {
        for (auto& source : _sources) source->get_parameters(params, filter);
        for (auto& psfmodel : _psfmodels) psfmodel->get_parameters(params, filter);
        _data->get_parameters(params, filter);
        return params;
    }

    ParamCRefs& get_parameters_const(ParamCRefs& params, ParamFilter* filter = nullptr) const override {
        for (auto& source : _sources) source->get_parameters_const(params, filter);
        for (auto& psfmodel : _psfmodels) psfmodel->get_parameters_const(params, filter);
        _data->get_parameters_const(params, filter);
        return params;
    }

    /// Same as get_parameters(), but for a single Observation with index idx in _data
    ParamRefs& get_parameters_observation(ParamRefs& params, size_t idx,
                                          ParamFilter* filter = nullptr) const {
        _check_obs_idx(idx);
        for (auto& source : _sources) source->get_parameters(params, filter);
        _psfmodels[idx]->get_parameters(params, filter);
        _data->at(idx).get().get_parameters(params, filter);
        return params;
    }

    /// Same as get_parameters_const(), but for a single Observation with index idx in _data
    ParamCRefs& get_parameters_observation_const(ParamCRefs& params, size_t idx,
                                                 ParamFilter* filter = nullptr) const {
        _check_obs_idx(idx);
        for (auto& source : _sources) source->get_parameters_const(params, filter);
        _psfmodels[idx]->get_parameters_const(params, filter);
        _data->at(idx).get().get_parameters_const(params, filter);
        return params;
    }

    /// Return _priors, the list of Prior instances
    Priors get_priors() const { return _priors; }

    /// Return _psfmodels, the list of PsfModel instances for each Observation in _data
    PsfModels get_psfmodels() const { return _psfmodels; }

    /// Return _sources, the list of Source instances for each Observation in _data
    Sources get_sources() const { return _sources; }

    void set_extra_param_factors(const Channel& channel, ExtraParamFactors& factors,
                                 size_t index) const override {
        for (auto& source : _sources) {
            source->set_extra_param_factors(channel, factors, index);
            index += source->get_n_gaussians(channel);
        }
    }

    void set_grad_param_factors(const Channel& channel, GradParamFactors& factors,
                                size_t index) const override {
        for (auto& source : _sources) {
            source->set_grad_param_factors(channel, factors, index);
            index += source->get_n_gaussians(channel);
        }
    }

    /**
     * Setup Evaluator instances for every Observation in _data using the given EvaluatorMode.
     *
     * @param mode The EvaluatorMode to use for all Evaluator instances.
     * @param outputs A vector of vectors of Image outputs for each Evaluator (created if empty and needed).
     * @param residuals A vector of residual Images for each Evaluator (created if empty and needed).
     * @param outputs_prior A vector of prior output (Jacobian) Images for each Evaluator
     *                      (created if empty and needed).
     * @param residuals_prior A vector of prior residual Images for each Evaluator
     *                             (created if empty and needed).
     * @param force Whether to force setting up even if already set up in the same mode
     * @param print Whether to print diagnostic statements to stdout.
     *
     * @note Different modes require different sized vectors for outputs
     *       EvaluatorMode::jacobian requires one Image per free parameter per Observation.
     *       EvaluatorMode::loglike_grad requires only one Image(n_rows=1, n_cols=n_params_free)
     */
    void setup_evaluators(EvaluatorMode mode = EvaluatorMode::image,
                          std::vector<std::vector<std::shared_ptr<Image>>> outputs = {},
                          std::vector<std::shared_ptr<Image>> residuals = {},
                          std::vector<std::shared_ptr<Image>> outputs_prior = {},
                          std::shared_ptr<Image> residuals_prior = nullptr, bool force = false,
                          bool print = false) {
        const size_t n_outputs = outputs.size();
        const bool has_outputs = n_outputs > 0;
        if (has_outputs) {
            if (n_outputs != size()) {
                throw std::invalid_argument("outputs.size()=" + std::to_string(n_outputs)
                                            + "!=this->size()=" + std::to_string(size()));
            }
        }
        const size_t n_residuals = residuals.size();
        const bool has_residuals = n_residuals > 0;
        if (has_residuals) {
            if (n_residuals != size()) {
                throw std::invalid_argument("residuals.size()=" + std::to_string(n_residuals)
                                            + "!=this->size()=" + std::to_string(size()));
            }
        }
        if (force || (!_is_setup || (mode != _mode))) {
            _evaluators.clear();
            _grads.clear();
            _outputs.clear();
            _outputs_prior.clear();
            _residuals_prior = nullptr;
            _evaluators.reserve(size());
            _is_setup = false;

            const bool is_jacobian = mode == EvaluatorMode::jacobian;
            const bool is_loglike_grad = mode == EvaluatorMode::loglike_grad;

            std::vector<ParamBaseCRef> params_free = {};
            if (is_jacobian || is_loglike_grad) {
                auto filter_free = g2f::ParamFilter{false, true, true, true};
                this->get_parameters_const(params_free, &filter_free);
                params_free = nonconsecutive_unique<ParamBaseCRef>(params_free);
                _n_params_free = params_free.size();
                if (print) {
                    std::cout << "n_params_free=" << _n_params_free << " from params_free:" << std::endl;
                    std::cout << str_iter_ref<true>(params_free) << std::endl;
                }

                _grads.reserve(size());
                _offsets_params.clear();
            }

            size_t n_residuals_prior = 0;
            bool do_residuals_prior = (_size_priors > 0)
                                      && (is_loglike_grad || is_jacobian || (mode == EvaluatorMode::loglike)
                                          || (_mode == EvaluatorMode::loglike_image));
            if (do_residuals_prior) {
                for (const auto& prior : this->_priors) {
                    n_residuals_prior += prior->size();
                }

                _residuals_prior
                        = residuals_prior == nullptr
                                  ? (_residuals_prior = std::make_shared<Image>(1, n_residuals_prior))
                                  : std::move(residuals_prior);
                if ((_residuals_prior->get_n_rows() != 1)
                    || (_residuals_prior->get_n_cols() != n_residuals_prior)) {
                    throw std::invalid_argument("residuals_prior=" + _residuals_prior->str()
                                                + " rows,cols != 1," + std::to_string(n_residuals_prior));
                }
                if (print) std::cout << "residuals_prior built/validated" << std::endl;
                for (size_t col = 0; col < n_residuals_prior; ++col) {
                    _residuals_prior->set_value(0, col, 0);
                }
                if (print) std::cout << "residuals_prior reset" << std::endl;

                if (is_jacobian) {
                    if (print) {
                        std::cout << "setup_evaluators jacobian called for model:" << std::endl;
                        std::cout << this->str() << std::endl;
                    }

                    const size_t n_params_jac = _n_params_free + 1;
                    if (outputs_prior.size() == 0) {
                        for (size_t idx_jac = 0; idx_jac < n_params_jac; ++idx_jac) {
                            _outputs_prior.emplace_back(std::make_shared<Image>(1, n_residuals_prior));
                        }
                        if (print) std::cout << "outputs_prior built" << std::endl;
                    } else if (outputs_prior.size() != n_params_jac) {
                        throw std::invalid_argument(
                                "jacobian outputs_prior->size()=" + std::to_string(outputs_prior.size())
                                + " != (n_params_free + 1 = " + std::to_string(n_params_jac) + ")");
                    }
                    _outputs_prior.reserve(n_params_jac);
                    size_t idx_jac = 0;
                    for (auto& output_prior : outputs_prior) {
                        if (output_prior == nullptr) {
                            throw std::invalid_argument("output_prior[" + std::to_string(idx_jac)
                                                        + "] is null");
                        }
                        if ((output_prior->get_n_rows() != 1)
                            || (output_prior->get_n_cols() != n_residuals_prior)) {
                            throw std::invalid_argument("outputs_prior[" + std::to_string(idx_jac)
                                                        + "]=" + output_prior->str() + " rows,cols != 1,"
                                                        + std::to_string(n_residuals_prior));
                        }
                        for (size_t col = 0; col < n_residuals_prior; ++col) {
                            output_prior->set_value(0, col, 0);
                        }
                        _outputs_prior.emplace_back(std::move(output_prior));
                        idx_jac++;
                    }
                    if (print) std::cout << "outputs_prior validated" << std::endl;
                }
            }
            if (!is_jacobian) {
                _outputs.reserve(size());
            }

            const auto& channels = _data->get_channels();

            if (print) std::cout << "making evaluators" << std::endl;

            for (size_t idx = 0; idx < _size; ++idx) {
                auto result = this->_make_evaluator(idx, mode, outputs,
                                                    has_residuals ? residuals[idx] : nullptr, print);
                _evaluators.emplace_back(std::move(result.first));
                if (is_jacobian || is_loglike_grad) {
                    _grads.emplace_back(std::move(result.second.second));
                } else {
                    _outputs.emplace_back(std::move(result.second.first));
                }
            }

            if (print) std::cout << "evaluators made" << std::endl;

            std::set<ParamBaseCRef> params_prior;
            size_t idx_prior = 0;

            // Check that the priors are reasonable.
            for (const auto& prior : _priors) {
                auto eval_prior = prior->evaluate(is_jacobian);
                size_t n_resid = prior->size();
                size_t n_resid_eval = eval_prior.residuals.size();
                if (n_resid != eval_prior.residuals.size()) {
                    throw std::logic_error(prior->str() + ".size()=" + std::to_string(n_resid)
                                           + " != evaluate(true).residuals.size()="
                                           + std::to_string(n_resid_eval));
                }
                if (is_jacobian) {
                    size_t idx_param = 0;
                    for (const auto& param_jacob : eval_prior.jacobians) {
                        const auto& param = param_jacob.first;
                        if (params_prior.find(param) != params_prior.end()) {
                            throw std::invalid_argument(
                                    "_priors[" + std::to_string(idx_prior) + "]=" + prior->str()
                                    + ".evaluate(true)[" + std::to_string(idx_param)
                                    + "]=" + param.get().str() + " already in prior param set."
                                    + " Did you apply multiple priors to the same parameter?");
                        }
                        idx_param++;
                    }
                }
                idx_prior++;
            }
            if (print && (_priors.size() > 0)) std::cout << "priors validated" << std::endl;
            if (_size_priors > 0) {
                if ((_residuals_prior != nullptr) && (_residuals_prior->get_n_cols() != n_residuals_prior)) {
                    throw std::invalid_argument("residuals_prior=" + _residuals_prior->str() + " n_residuals="
                                                + std::to_string(n_residuals_prior) + "!= get_n_cols()="
                                                + std::to_string(_residuals_prior->get_n_cols()));
                }
            }
            _is_setup = true;
            _mode = mode;
        }
    }

    /// Get the size of this->_data
    size_t size() const { return _size; }

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override {
        std::string str = type_name_str<Model>(false, namespace_separator) + "("
                          + (name_keywords ? "data=" : "") + _data->repr(name_keywords, namespace_separator)
                          + ", " + (name_keywords ? "psfmodels=[" : "[");
        for (const auto& x : _psfmodels) str += x->repr(name_keywords, namespace_separator) + ",";
        str += (name_keywords ? "], sources=[" : "], [");
        for (const auto& s : _sources) str += s->repr(name_keywords, namespace_separator) + ",";
        str += (name_keywords ? "], priors=[" : "], [");
        return str + "])";
    }

    std::string str() const override {
        std::string str = type_name_str<Model>(true) + "(data=" + _data->str() + ", psfmodels=[";
        for (const auto& x : _psfmodels) str += x->str() + ",";
        str += "], sources=[";
        for (const auto& x : _sources) str += x->str() + ",";
        str += "], priors=[";
        for (const auto& x : _priors) str += x->str() + ",";
        return str + "])";
    }

    /**
     * Verify that the Jacobian is correct by comparing to finite differences.
     *
     * @param findiff_frac The value of the finite difference increment, as a fraction of the parameter value.
     * @param findiff_add The minimum value of the finite difference increment.
     * @param rtol The allowed relative tolerance in the Jacobian as compared to the finite difference.
     * @param atol The allowed absolute tolerance in the Jacobian as compared to the finite difference.
     * @return A vector of error messages, one for each Parameter that failed verification.
     *
     * @note Verification is done using isclose(), which is modelled after Python's numpy.isclose.
     */
    std::vector<std::string> verify_jacobian(double findiff_frac = 1e-5, double findiff_add = 1e-5,
                                             double rtol = 1e-3, double atol = 1e-3) {
        if (_mode != EvaluatorMode::jacobian) {
            this->setup_evaluators(EvaluatorMode::jacobian);
        }

        const size_t n_obs = this->size();
        ParamFilter filter{false, true};
        std::vector<std::string> errors;

        for (size_t idx = 0; idx < n_obs; ++idx) {
            const auto& grads = *(_grads.at(idx));
            const auto& observation = _data->at(idx).get();
            const auto& sigma_inv = observation.get_sigma_inverse();

            const Channel& channel = observation.get_channel();
            const size_t n_cols = observation.get_n_cols();
            const size_t n_rows = observation.get_n_rows();

            filter.channel = channel;

            ParamRefs params{};
            this->get_parameters_observation(params, idx, &filter);
            params = nonconsecutive_unique(params);

            size_t idx_param_max = 0;
            for (const auto& param : params) {
                auto found = _offsets_params.find(param);
                if (found == _offsets_params.end()) {
                    throw std::runtime_error(
                            param.get().str()
                            + " not found in offsets; was it freed after setup_evaluators was called?");
                }
                size_t idx_param = found->second;
                if (idx_param > idx_param_max) idx_param_max = idx_param;
            }

            if (!(idx_param_max < grads.size())) {
                throw std::runtime_error("idx_param_max=" + std::to_string(idx_param_max)
                                         + " !< grads.size()=" + std::to_string(grads.size())
                                         + "; is the jacobian array large enough?");
            }

            double loglike_jac = this->_evaluate_observation(idx);
            const auto& grad = grads[0];
            size_t n_nonzero = 0;
            for (unsigned int i = 0; i < n_cols; ++i) {
                for (unsigned int j = 0; j < n_rows; ++j) {
                    n_nonzero += grad.get_value(j, i) != 0;
                }
            }
            if (n_nonzero > 0) {
                errors.push_back("n_nonzero grads[0] entries=" + std::to_string(n_nonzero));
            }

            auto result = _make_evaluator(idx, EvaluatorMode::loglike_image);
            double loglike_img = result.first->loglike_pixel();
            if (loglike_jac != loglike_img) {
                errors.push_back("loglike_jac=" + to_string_float(loglike_jac) + "loglike_img="
                                 + to_string_float(loglike_img) + "for obs[" + std::to_string(idx) + "]:");
            }

            const auto& image_current = *(result.second.first);

            for (auto& param_ref : params) {
                auto& param = param_ref.get();
                const size_t idx_param = _offsets_params.find(param_ref)->second;
                const auto& grad = grads[idx_param];

                const double value_init = param.get_value();
                const double value_transformed = param.get_value_transformed();
                double diff = value_transformed * findiff_frac;
                if (std::abs(diff) < findiff_add) diff = findiff_add;
                diff = finite_difference_param(param, diff);

                auto result_diff = _make_evaluator(idx, EvaluatorMode::image);

                result_diff.first->loglike_pixel();

                size_t n_failed = 0;
                std::vector<double> ratios = {};
                for (unsigned int i = 0; i < n_cols; ++i) {
                    for (unsigned int j = 0; j < n_rows; ++j) {
                        double delta = sigma_inv.get_value(j, i) / diff
                                       * (result_diff.second.first->get_value(j, i)
                                          - image_current.get_value(j, i));
                        double jac = grad.get_value(j, i);

                        auto close = isclose(jac, delta, rtol, atol);
                        if (!close.isclose) {
                            n_failed++;
                            ratios.push_back(jac / delta);
                        }
                    }
                }

                param.set_value(value_init);
                const double value_new = param.get_value();
                if (value_new != value_init) {
                    throw std::logic_error("Could not return param=" + param.str()
                                           + " to value_init=" + to_string_float(value_init) + "; diff="
                                           + to_string_float(value_new - value_init) + "); check limits");
                }
                if (n_failed > 0) {
                    std::sort(ratios.begin(), ratios.end());
                    double median = ratios[ratios.size() / 2];
                    errors.push_back("Param[" + std::to_string(idx_param) + "]=" + param.str()
                                     + " failed for obs[" + std::to_string(idx)
                                     + "]: n_failed=" + std::to_string(n_failed)
                                     + "; median evaluated/expected=" + std::to_string(median));
                }
            }
        }
        filter.channel = std::nullopt;

        ParamRefs params{};
        this->get_parameters_observation(params, 0, &filter);
        params = nonconsecutive_unique(params);

        for (const auto& prior : this->_priors) {
            auto result_base = prior->evaluate(true);
            auto residuals_base = result_base.residuals;
            const size_t n_values = residuals_base.size();

            for (const auto& [param_cref, values_base] : result_base.jacobians) {
                if (values_base.size() != n_values) {
                    throw std::logic_error("values_base.size()=" + std::to_string(values_base.size())
                                           + " != residuals_base.size()=" + std::to_string(n_values));
                }

                // TODO: Determine if it's possible to get mutable ref from set.find
                // (set.find only ever seems to return cref, so maybe not)
                auto found = std::find(params.begin(), params.end(), param_cref);
                if (found == params.end()) {
                    throw std::logic_error("Couldn't find " + param_cref.get().str() + " in all params");
                }
                auto& param = (*found).get();

                const double value_init = param.get_value();
                const double value_transformed = param.get_value_transformed();
                double delta = value_transformed * findiff_frac;
                if (std::abs(delta) < findiff_add) delta = findiff_add;
                delta = finite_difference_param(param, delta);
                auto result_new = prior->evaluate(true);

                const size_t idx_param = _offsets_params.at(param);

                std::vector<double> ratios = {};
                size_t n_failed = 0;
                size_t idx_value = 0;
                for (const double jac_base : values_base) {
                    double jac_findiff
                            = (result_new.residuals[idx_value] - result_base.residuals[idx_value]) / delta;
                    auto close = isclose(jac_base, jac_findiff, rtol, atol);
                    if (!close.isclose) {
                        n_failed++;
                        prior->evaluate(true);
                        ratios.push_back(jac_base / jac_findiff);
                    }
                    idx_value++;
                }
                param.set_value(value_init);
                const double value_new = param.get_value();
                if (value_new != value_init) {
                    throw std::logic_error("Could not return param=" + param.str()
                                           + " to value_init=" + to_string_float(value_init) + "; diff="
                                           + to_string_float(value_new - value_init) + "); check limits");
                }
                if (n_failed > 0) {
                    std::sort(ratios.begin(), ratios.end());
                    double median = ratios[ratios.size() / 2];
                    errors.push_back("Param[" + std::to_string(idx_param) + "]=" + param.str()
                                     + " failed for prior=" + prior->str()
                                     + ": n_failed=" + std::to_string(n_failed)
                                     + "; median evaluated/expected=" + std::to_string(median));
                }
            }
        }
        return errors;
    }
};

}  // namespace lsst::gauss2d::fit

#endif
