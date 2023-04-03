#ifndef GAUSS2D_FIT_MODEL_H
#define GAUSS2D_FIT_MODEL_H

#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "gauss2d/evaluate.h"

#include "data.h"
#include "gauss2d/image.h"
#include "parametricmodel.h"
#include "param_filter.h"
#include "psfmodel.h"
#include "source.h"
#include "util.h"

namespace gauss2d::fit
{

//enum class Renderer { gauss2d };

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
class Model : public ParametricModel
{
public:
    typedef GaussianEvaluator<T, Image, Indices> Evaluator;
    typedef std::map<std::reference_wrapper<const Channel>,
        std::vector<std::unique_ptr<const Gaussians>>> GaussiansMap;
    typedef Data<T, Image, Mask> ModelData;
    typedef typename ModelData::Observation Observation;
    typedef std::vector<std::shared_ptr<PsfModel>> PsfModels;
    typedef std::vector<std::shared_ptr<Source>> Sources;

    /// Valid forms of evaluation to setup GaussianEvaluator instances for.
    enum class EvaluatorMode {
        image,         ///< Compute the model images only
        loglike,       ///< Compute the log likelihood only
        loglike_image, ///< Compute the model images and log likelihood
        loglike_grad,  ///< Compute the gradients of the log likelihood
        jacobian       ///< Compute the model Jacobian
    };

private:
    std::shared_ptr<const ModelData> _data;
    std::vector<std::unique_ptr<Evaluator>> _evaluators = {};
    ExtraParamFactors _factors_extra = {};
    std::vector<std::weak_ptr<Image>> _factors_extra_in;
    GradParamFactors _factors_grad = {};
    std::vector<std::weak_ptr<Image>> _factors_grad_in;
    GaussiansMap _gaussians_srcs;
    std::vector<std::shared_ptr<ImageArray<double, Image>>> _grads;
    bool _is_setup = false;
    std::vector<std::weak_ptr<Indices>> _map_extra_in;
    std::vector<std::weak_ptr<Indices>> _map_grad_in;
    EvaluatorMode _mode;
    std::vector<std::shared_ptr<Image>> _outputs = {};
    std::vector<std::unique_ptr<const gauss2d::Gaussians>> _psfcomps;
    PsfModels _psfmodels;
    size_t _size;
    Sources _sources;

    void _check_obs_idx(size_t idx) const {
        const size_t n_obs = _data->size();
        if(!(idx < n_obs)) throw std::invalid_argument("idx=" + std::to_string(idx) + " !< data->size()="
            + std::to_string(n_obs));
    }

    /// Return a map of source Gaussians for each Channel in Data
    GaussiansMap _get_gaussians_srcs() const {
        GaussiansMap gaussians_srcs{};
        size_t n_gaussians_src = 0;

        for(const Channel & channel : _data->get_channels()) {
            size_t n_gaussians_c = 0;
            gaussians_srcs[channel] = std::vector<std::unique_ptr<const Gaussians>>();
            auto & gaussians_src = gaussians_srcs[channel];
            for(const auto & src : this->get_sources()) {
                gaussians_src.push_back(src->get_gaussians(channel));
                n_gaussians_c += gaussians_src.back()->size();
            }
            if(n_gaussians_src == 0) {
                n_gaussians_src = n_gaussians_c;
            } else if(n_gaussians_src != n_gaussians_c) {
                throw std::logic_error("n_gaussians[" + channel.str() + "]=" + std::to_string(n_gaussians_c)
                    + "!=n_gaussians=" + std::to_string(n_gaussians_src));
            }
        }
        return gaussians_srcs;
    }

    /// Evaluate a single observation using the current EvaluatorMode in _mode
    double _evaluate_observation(size_t idx) {
        _check_obs_idx(idx);
        if(_mode == EvaluatorMode::jacobian) {
            const auto & channel = _data->at(idx).get().get_channel();
            const size_t n_gaussians_psf = _psfcomps[idx]->size();
            const auto & gaussians_src = _gaussians_srcs.at(channel);

            const auto & factors_extra_in = _factors_extra_in[idx].lock();
            const auto & factors_grad_in = _factors_grad_in[idx].lock();

            size_t offset = 0;
            for(size_t i_psf = 0; i_psf < n_gaussians_psf; ++i_psf)
            {
                size_t i_src = 0;
                for(const auto & src : this->get_sources())
                {
                    const size_t n_gauss_src = gaussians_src[i_src++]->size();
                    if(n_gauss_src != src->get_n_gaussians(channel)) {
                        throw std::logic_error(this->str() + "._data[" + channel.str()
                            + "].get_n_gaussians(channel)=" + std::to_string(src->get_n_gaussians(channel))
                            + " != get_gaussians(channel).size()=" + std::to_string(n_gauss_src));
                    }
                    src->set_grad_param_factors(channel, _factors_grad, offset);
                    src->set_extra_param_factors(channel, _factors_extra, offset);
                    // Copy values to the actual input arrays used by the evaluator
                    const size_t row_max = offset + n_gauss_src;
                    for(size_t row = offset; row < row_max; ++row) {
                        for(size_t col = 0; col < gauss2d::N_EXTRA_FACTOR; ++col) {
                            factors_extra_in->set_value_unchecked(row, col, _factors_extra[row][col]);
                        }
                        for(size_t col = 0; col < gauss2d::N_PARAMS_GAUSS2D; ++col) {
                            factors_grad_in->set_value_unchecked(row, col, _factors_grad[row][col]);
                        }
                    }
                    offset = row_max;
                }
            }
        }
        return _evaluators[idx]->loglike_pixel();
    }

    /**
     * Setup GaussianEvaluator instances for every Observation in Data using one EvaluatorMode.
     *
     * @param idx_obs The number index of the Observation in Data.
     * @param mode The EvaluatorMode for Evaluator instances.
     * @param outputs The list of Image instances to output to, ordered as Data.
     * @param residual The residual Image for the Evaluator to use.
     * @param print Whether to print diagnostic information to stdout.
     * @return A pair of Evaluator instances, along with the associated output Image
     *         and gradient ImageArray.
     *
     * @note Different modes require different lengths of outputs.
     * @note Not all EvaluatorMode options are currently implemented.
     */
    std::pair<
        std::unique_ptr<Evaluator>,
        std::pair<std::shared_ptr<Image>, std::shared_ptr<ImageArray<double, Image>>>
    > _make_evaluator(
        const size_t idx_obs,
        EvaluatorMode mode=EvaluatorMode::image,
        std::vector<std::shared_ptr<Image>> outputs={},
        std::shared_ptr<Image> residual=nullptr,
        bool print=false
    ) {
        _check_obs_idx(idx_obs);
        // TODO: implement
        if(mode == EvaluatorMode::loglike_grad) throw std::runtime_error(
            "loglike_grad mode not implemented yet");
        const Observation & observation = _data->at(idx_obs).get();
        const Channel & channel = observation.get_channel();
        const auto & psfcomps = _psfcomps[idx_obs];
        const size_t n_c = observation.get_n_cols();
        const size_t n_r = observation.get_n_rows();

        const size_t n_outputs = outputs.size();
        const bool has_outputs = n_outputs > 0;
        if(has_outputs) {
            // TODO: Enable for image and loglike_image
            if(mode != EvaluatorMode::jacobian) {
                throw std::invalid_argument("outputs not implemented for non-jacobian modes");
            }
        }

        const bool is_mode_image = mode == EvaluatorMode::image;
        const bool is_mode_jacobian = mode == EvaluatorMode::jacobian;
        const bool do_output = is_mode_image || (mode == EvaluatorMode::loglike_image);

        std::shared_ptr<const CoordinateSystem> coordsys = observation.get_image().get_coordsys_ptr_const();
        std::shared_ptr<Image> output = !do_output ? nullptr : 
            std::make_shared<Image>(n_r, n_c, coordsys);

        auto data_eval = !is_mode_image ? observation.get_image_ptr_const() : nullptr;
        auto sigma_inv = !is_mode_image ? observation.get_sigma_inverse_ptr_const() : nullptr;

        const auto & gaussians_src = _gaussians_srcs.at(channel);
        const size_t n_gaussians_psf = psfcomps->size();

        for(size_t p = 0; p < n_gaussians_psf; ++p)
        {
            const auto & psfcomp = psfcomps->at(p);
            const auto value = psfcomp.get_integral_value();
            // TODO: Downgrade to warning and/or handle neatly
            if(!(value > 0))
            {
                throw std::runtime_error("psfcomps[" + std::to_string(p) + "]=" + psfcomp.str()
                    + " get_integral_value()=" + std::to_string(value) + "!>0");
            }
        }

        size_t n_gaussians_src = 0;
        for(const auto & gaussians : gaussians_src) n_gaussians_src += gaussians->size();
        const size_t n_gaussians_conv = n_gaussians_src*n_gaussians_psf;

        ConvolvedGaussians::Data data {};

        for(auto psfit = psfcomps->cbegin(); psfit < psfcomps->cend(); psfit++) {
            for(const auto & gaussians : gaussians_src) {
                for(const auto & gaussian : *gaussians) {
                    data.emplace_back(std::make_shared<ConvolvedGaussian>(gaussian, *psfit));
                }
            }
        }
        auto gaussians = std::make_shared<const ConvolvedGaussians>(data);
        if(gaussians->size() != n_gaussians_conv) {
            throw std::logic_error("gaussians.size()=" + std::to_string(gaussians->size())
                + "!= n_gaussians_conv=" + std::to_string(n_gaussians_conv) + "= n_gaussians_src="
                + std::to_string(n_gaussians_src) + "= n_gaussians_psf=" + std::to_string(n_gaussians_psf)
                + "(is_mode_jacobian=" + std::to_string(is_mode_jacobian)
            );
        }

        std::shared_ptr<const Indices> map_extra_in, map_grad_in;
        std::shared_ptr<const Image> factors_extra_in, factors_grad_in;

        if(is_mode_jacobian) {
            std::weak_ptr<Image> factors_extra_weak = _factors_extra_in[idx_obs];
            std::weak_ptr<Image> factors_grad_weak = _factors_grad_in[idx_obs];
            std::weak_ptr<Indices> map_extra_weak = _map_extra_in[idx_obs];
            std::weak_ptr<Indices> map_grad_weak = _map_grad_in[idx_obs];

            if(print) {
                this->_make_param_maps<true>(
                    factors_extra_weak, factors_grad_weak, map_extra_weak, map_grad_weak,
                    n_gaussians_conv, n_gaussians_psf, channel, coordsys, idx_obs,
                    map_extra_in, map_grad_in, factors_extra_in, factors_grad_in
                );
            } else {
                this->_make_param_maps<false>(
                    factors_extra_weak, factors_grad_weak, map_extra_weak, map_grad_weak,
                    n_gaussians_conv, n_gaussians_psf, channel, coordsys, idx_obs,
                    map_extra_in, map_grad_in, factors_extra_in, factors_grad_in
                );
            }
        }
        
        std::shared_ptr<ImageArray<double, Image>> grads;
        if(is_mode_jacobian) {
            typename ImageArray<double, Image>::Data arrays {};
            ParamCRefs params_free {};
            auto filter_free = g2f::ParamFilter{false, true, true, true, channel};
            this->get_parameters_observation_const(params_free, idx_obs, &filter_free);
            params_free = nonconsecutive_unique<ParamBaseCRef>(params_free);

            const size_t n_free = params_free.size() + 1;

            if(has_outputs) {
                const size_t size_jac = outputs.size();
                if(size_jac != n_free) {
                    throw std::invalid_argument(
                        "outputs[" + std::to_string(idx_obs) + "].size()=" + std::to_string(size_jac)
                            + "!n_free=" + std::to_string(n_free)
                    );
                }
            }

            for(size_t idx_img = 0; idx_img < n_free; ++idx_img) {
                std::shared_ptr<Image> image;
                if(has_outputs) {
                    image = outputs[idx_img];
                    if(image == nullptr) {
                        throw std::invalid_argument("output[" + std::to_string(idx_img)
                            + "] can't be null for obs #" + std::to_string(idx_obs));
                    }
                    const auto & image_ref = observation.get_image();
                    if(!images_compatible(*image, image_ref)) {
                        throw std::invalid_argument("outputs[" + std::to_string(idx_img) + "]="
                            + image->str() + " incompatible with corresponding image=" + image_ref.str()
                            + "for obs #" + std::to_string(idx_obs));
                    }
                } else {
                    image = std::make_shared<Image>(n_r, n_c, coordsys);
                }
                arrays.emplace_back(image);
            }
            grads = std::make_shared<ImageArray<double, Image>>(&arrays);
        } else {
            grads = nullptr;
        }
        if(print) {
            std::cout << observation.str() << std::endl;
            std::cout << gaussians->str() << std::endl;
        }

        auto grads2 = grads;
        auto output2 = output;
	    std::pair<std::shared_ptr<Image>, std::shared_ptr<ImageArray<double, Image>>> second = {
            std::move(output2), std::move(grads2)
        };

        return {
            std::make_unique<Evaluator>(
                gaussians,
                coordsys,
                data_eval,
                sigma_inv,
                output,
                residual,
                grads,
                map_grad_in,
                factors_grad_in,
                map_extra_in,
                factors_extra_in
                //, background
            ),
            std::move(second)
        };
    }

    /// Make new factor/map inputs for GaussianEvaluator, if needed.
    template <bool print>
    void _make_param_maps(
        std::weak_ptr<Image> factors_extra_weak,
        std::weak_ptr<Image> factors_grad_weak,
        std::weak_ptr<Indices> map_extra_weak,
        std::weak_ptr<Indices> map_grad_weak,
        size_t n_gaussians_conv,
        size_t n_gaussians_psf,
        const Channel & channel,
        std::shared_ptr<const CoordinateSystem> coordsys,
        size_t idx_obs,
        std::shared_ptr<const Indices> & map_extra_in,
        std::shared_ptr<const Indices> & map_grad_in,
        std::shared_ptr<const Image> & factors_extra_in,
        std::shared_ptr<const Image> & factors_grad_in
    ) {
        unsigned int n_expired = map_extra_weak.expired() + map_grad_weak.expired() + 
        factors_grad_weak.expired() + factors_extra_weak.expired();
        const bool expired = n_expired == 4;
        if(!((n_expired == 0) || expired)) {
            throw std::logic_error("jacobian n_expired=" + std::to_string(n_expired) + " not in (0,4)");
        }
        auto map_extra_mut = expired ? std::make_shared<Indices>(n_gaussians_conv, 2, coordsys)
            : map_extra_weak.lock();
        auto map_grad_mut = expired ? std::make_shared<Indices>(
            n_gaussians_conv, gauss2d::N_PARAMS_GAUSS2D, coordsys) : map_grad_weak.lock();
        auto factors_extra_mut = expired ? std::make_shared<Image>(n_gaussians_conv, 3, coordsys)
            : factors_extra_weak.lock();
        auto factors_grad_mut = expired ? std::make_shared<Image>(
            n_gaussians_conv, gauss2d::N_PARAMS_GAUSS2D, coordsys) : factors_grad_weak.lock();;

        ExtraParamMap map_extra = {};
        GradParamMap map_grad = {};
        ExtraParamFactors & factors_extra = _factors_extra;
        GradParamFactors & factors_grad = _factors_grad;
        factors_extra.resize(0);
        factors_grad.resize(0);

        map_extra.reserve(n_gaussians_conv);
        map_grad.reserve(n_gaussians_conv);
        factors_extra.reserve(n_gaussians_conv);
        factors_grad.reserve(n_gaussians_conv);

        ParameterMap offsets {};
        /*
            Make (and set) the param map; only make the factors (re-evaluated each time) 
            The maps depend on the order of free parameters.
            the factors depend on their values, and can change every evaluation.
        */
        for(size_t i_psf = 0; i_psf < n_gaussians_psf; ++i_psf)
        {
            for(const auto & src : this->get_sources())
            {
                src->add_grad_param_map(channel, map_grad, offsets);
                src->add_extra_param_map(channel, map_extra, map_grad, offsets);

                src->add_grad_param_factors(channel, factors_grad);
                src->add_extra_param_factors(channel, factors_extra);
            }
        }
        if(map_grad.size() != n_gaussians_conv) {
            throw std::logic_error(
                "map_grad.size()=" + std::to_string(map_grad.size()) + "!= n_gaussians_conv="
                + std::to_string(n_gaussians_conv)
            );
        }
        /*
            Validate the map if there is already another Jacobian evaluator.

            Do not support multiple evaluators set up with different maps (yet).
        */
        if(!expired) {
            std::string sizes_mut = std::to_string(map_extra_mut->size()) + ","
                + std::to_string(map_grad_mut->size()) + ","
                + std::to_string(factors_extra_mut->size()) + ","
                + std::to_string(factors_grad_mut->size());

            std::string sizes = std::to_string(map_extra.size()) + ","
                + std::to_string(map_grad.size()) + ","
                + std::to_string(factors_extra.size()) + ","
                + std::to_string(factors_grad.size());
            
            if(sizes != sizes_mut) throw std::runtime_error(
                "make_evaluator Jacobian map_extra,map_grad,factors_extra,factors_grad sizes_new="
                + sizes + "!=sizes_mut=" + sizes_mut + "; did you make a new evaluator with different"
                "free parameters from an old one?");
            
            std::array<std::string, 4> errmsgs = {"", "", "", ""};

            size_t idx_err = 0;
            if constexpr (print) std::cout << "map_extra=[" << std::endl;
            for(size_t idx_row = 0; idx_row < map_extra.size(); idx_row++) {
                const auto & row = map_extra[idx_row];
                for(size_t col=0; col < gauss2d::N_EXTRA_MAP; col++) {
                    const auto value_mut = map_extra_mut->get_value_unchecked(idx_row, col);
                    if((idx_err < 2) && (row[col] != value_mut)) errmsgs[idx_err++] = 
                        "map_extra[" + std::to_string(idx_row) + "][" + std::to_string(col)
                        + "]=" + std::to_string(row[col]) + " != old="
                        + std::to_string(value_mut) + "\n";
                }
                if constexpr (print) {
                    std::cout << "    ";
                    stream_iter_ref(row, std::cout);
                    std::cout << std::endl;
                }
            }
            if constexpr (print) std::cout << "]" << std::endl << "map_grad=[" << std::endl;
            for(size_t idx_row = 0; idx_row < map_grad.size(); idx_row++) {
                const auto & row = map_grad[idx_row];
                for(size_t col=0; col < gauss2d::N_EXTRA_MAP; col++) {
                    const auto value_mut = map_grad_mut->get_value_unchecked(idx_row, col);
                    if((idx_err < 4) && (row[col] != value_mut)) errmsgs[idx_err++] = 
                        "map_grad[" + std::to_string(idx_row) + "][" + std::to_string(col)
                        + "]=" + std::to_string(row[col]) + " != old="
                        + std::to_string(value_mut) + "\n";
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
            if(!errmsg_extra.empty() || !errmsg_grad.empty()) {
                std::string errmsg = "make_evaluator Jacobian map_extra/map_grad have different values"
                    " from old; did you make a new evaluator with different free parameters from an "
                    " old one? Errors:\n" + errmsg_extra;
                if(!errmsg_extra.empty() && !errmsg_grad.empty()) errmsg += "...\n";
                errmsg += errmsg_grad;
                throw std::runtime_error(errmsg);
            }
        } else {
            /*
                Maps are stored as vectors in gauss2dfit (to avoid unnecessary templating)
                The inputs to gauss2d's Evaluator are templated Images, so copy the values
            */
            if constexpr (print) std::cout << "factors/maps: {" << std::endl;
            for(size_t row = 0; row < n_gaussians_conv; ++row)
            {
                for(size_t col = 0; col < gauss2d::N_EXTRA_MAP; ++col) {
                    map_extra_mut->set_value_unchecked(row, col, (map_extra)[row][col]);
                }
                for(size_t col = 0; col < gauss2d::N_EXTRA_FACTOR; ++col) {
                    factors_extra_mut->set_value_unchecked(row, col, (factors_extra)[row][col]);
                }
                for(size_t col = 0; col < gauss2d::N_PARAMS_GAUSS2D; ++col) {
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
                    std::cout << std::endl;
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
    void add_extra_param_map(
            const Channel & channel, ExtraParamMap & map_extra, const GradParamMap & map_grad,
            ParameterMap & offsets
    ) const override
    {
        for(auto & source : _sources) source->add_extra_param_map(channel, map_extra, map_grad, offsets);
    }
    void add_extra_param_factors(const Channel & channel, ExtraParamFactors & factors) const override
    {
        for(auto & source : _sources) source->add_extra_param_factors(channel, factors);
    }
    void add_grad_param_map(const Channel & channel, GradParamMap & map, ParameterMap & offsets
        ) const override
    {
        for(auto & source : _sources) source->add_grad_param_map(channel, map, offsets);
    }
    void add_grad_param_factors(const Channel & channel, GradParamFactors & factors) const override
    {
        for(auto & source : _sources) source->add_grad_param_factors(channel, factors);
    }

    /// Evaluate the model for every Observation in _data.
    std::vector<double> evaluate() {
        if(!_is_setup) throw std::runtime_error("Can't call evaluate before setup_evaluators");

        std::vector<double> result(_evaluators.size());
        
        for(size_t idx = 0; idx < _size; ++idx) {
            result[idx] = _evaluate_observation(idx);
        }

        return result;
    }

    /// Evaluate a single observation with the given index in _data.
    double evaluate_observation(size_t idx) {
        _check_obs_idx(idx);
        return _evaluate_observation(_evaluators[idx]);
    }

    /// Return _data
    std::shared_ptr<const ModelData> get_data() const { return _data; }

    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override {
        std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
        // TODO: Rethink this; it's not sufficient unless sources are single component gaussians
        in.reserve(_sources.size());

        for(auto & source : _sources) {
            auto data = source->get_gaussians(channel)->get_data();
            in.emplace_back(data);
        }

        return std::make_unique<gauss2d::Gaussians>(in);
    }

    size_t get_n_gaussians(const Channel & channel) const override {
        size_t n_g = 0;
        for(auto & source : _sources) n_g += source->get_n_gaussians(channel);
        return n_g;
    }

    /// Return _outputs (output Image instances for each Observation in _data)
    std::vector<std::shared_ptr<Image>> get_outputs() const { return _outputs; }

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override
    {
        for(auto & source : _sources) source->get_parameters(params, filter);
        for(auto & psfmodel : _psfmodels) psfmodel->get_parameters(params, filter);
        _data->get_parameters(params, filter);
        return params;
    }

    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override
    {
        for(auto & source : _sources) source->get_parameters_const(params, filter);
        for(auto & psfmodel : _psfmodels) psfmodel->get_parameters_const(params, filter);
        _data->get_parameters_const(params, filter);
        return params;
    }

    /// Same as get_parameters(), but for a single Observation with index idx in _data
    ParamRefs & get_parameters_observation(ParamRefs & params, size_t idx, ParamFilter * filter = nullptr)
    const {
        _check_obs_idx(idx);
        for(auto & source : _sources) source->get_parameters(params, filter);
        _psfmodels[idx]->get_parameters(params, filter);
        _data->at(idx).get().get_parameters(params, filter);
        return params;
    }

    /// Same as get_parameters_const(), but for a single Observation with index idx in _data
    ParamCRefs & get_parameters_observation_const(ParamCRefs & params, size_t idx,
        ParamFilter * filter = nullptr) const {
        _check_obs_idx(idx);
        for(auto & source : _sources) source->get_parameters_const(params, filter);
        _psfmodels[idx]->get_parameters_const(params, filter);
        _data->at(idx).get().get_parameters_const(params, filter);
        return params;
    }

    /// Return _psfmodels, the list of PsfModel instances for each Observation in _data
    PsfModels get_psfmodels() const {
        return _psfmodels;
    }

    /// Return _sources, the list of Source instances for each Observation in _data
    Sources get_sources() const {
        return _sources;
    }

    void set_extra_param_factors(const Channel & channel, ExtraParamFactors & factors, size_t index)
        const override
    {
        for(auto & source : _sources) {
            source->set_extra_param_factors(channel, factors, index);
            index += source->get_n_gaussians(channel);
        }
    }

    void set_grad_param_factors(const Channel & channel, GradParamFactors & factors, size_t index)
        const override
    {
        for(auto & source : _sources) {
            source->set_grad_param_factors(channel, factors, index);
            index += source->get_n_gaussians(channel);
        }
    }

    /**
     * Setup Evaluator instances for every Observation in _data using the given EvaluatorMode.
     *
     * @param mode The EvaluatorMode to use for all Evaluator instances.
     * @param outputs A vector of vectors of Image outputs for each Evaluator (created if empty and needed).
     * @param residuals An array of residual
     * @param print Whether to print diagnostic statement to stdout.
     *
     * @note Different modes require different sized vectors for outputs
     *       EvaluatorMode::jacobian requires one Image per free parameter per Observation.
     */
    void setup_evaluators(
        EvaluatorMode mode=EvaluatorMode::image,
        std::vector<std::vector<std::shared_ptr<Image>>> outputs={},
        std::vector<std::shared_ptr<Image>> residuals={},
        bool print=false
    ) {
        const size_t n_outputs = outputs.size();
        const bool has_outputs = n_outputs > 0;
        if(has_outputs) {
            if(n_outputs != size()) {
                throw std::invalid_argument(
                    "outputs.size()=" + std::to_string(n_outputs)
                     + "!=this->size()=" + std::to_string(size())
                );
            }
        }
        const size_t n_residuals = residuals.size();
        const bool has_residuals = n_residuals > 0;
        if(has_residuals) {
            if(n_residuals != size()) {
                throw std::invalid_argument(
                    "residuals.size()=" + std::to_string(n_residuals)
                     + "!=this->size()=" + std::to_string(size())
                );
            }
        }
        if(!_is_setup || (mode != _mode)) {
            _evaluators.resize(0);
            _grads.resize(0);
            _outputs.resize(0);
            _evaluators.reserve(size());
            if(mode == EvaluatorMode::jacobian) {
                _grads.reserve(size());
            } else if(mode == EvaluatorMode::loglike_grad) {
                // TODO: implement
                throw std::runtime_error("loglike_grad mode not implemented yet");
            } else {
                _outputs.reserve(size());
            }
            
            const auto & channels = _data->get_channels();

            // Apparently this can't go directly in a ternary statement, so here it is
            std::vector<std::shared_ptr<Image>> outputs_obs{};

            for(size_t idx = 0; idx < _size; ++idx) {
                auto result = this->_make_evaluator(
                    idx, mode,
                    has_outputs ? outputs[idx] : outputs_obs,
                    has_residuals ? residuals[idx] : nullptr,
                    print
                );
                _evaluators.emplace_back(std::move(result.first));
                if(mode == EvaluatorMode::jacobian) {
                    _grads.emplace_back(std::move(result.second.second));
                } else if(mode == EvaluatorMode::loglike_grad) {
                    // TODO: implement
                    throw std::runtime_error("loglike_grad mode not implemented yet");
                } else {
                    _outputs.emplace_back(std::move(result.second.first));
                }
            }
            _is_setup = true;
            _mode = mode;
        }
    }

    /// Get the size of this->_data
    size_t size() const { return _size; }

    std::string repr(bool name_keywords = false) const override {
        std::string str = std::string("Model(") + (name_keywords ? "sources=[" : "[");
        for(const auto & s : _sources) str += s->repr(name_keywords) + ",";
        return str + "], " + (name_keywords ? "data=" : "") + _data->repr() + ")";
    }

    std::string str() const override {
        std::string str = "Model(sources=[";
        for(const auto & s : _sources) str += s->str() + ",";
        return str + "], data=" + _data->str() + ")";
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
    std::vector<std::string> verify_jacobian(
        double findiff_frac=1e-4,
        double findiff_add=1e-4,
        double rtol=1e-3,
        double atol=1e-3
    ) {
        if(_mode != EvaluatorMode::jacobian) {
            this->setup_evaluators(EvaluatorMode::jacobian);
            this->evaluate();
        }

        const size_t n_exp = _size;
        ParamFilter filter{false, true};
        std::vector<std::string> errors;

        for(size_t idx = 0; idx < n_exp; ++idx)
        {
            auto & evaluator = *(_evaluators.at(idx));
            const auto & grads = *(_grads.at(idx));
            const auto & observation = _data->at(idx).get();
            const auto & sigma_inv = observation.get_sigma_inverse();

            const Channel & channel = observation.get_channel();
            const size_t n_cols = observation.get_n_cols();
            const size_t n_rows = observation.get_n_rows();

            filter.channel = channel;

            ParamRefs params{};
            this->get_parameters_observation(params, idx, &filter);
            params = nonconsecutive_unique(params);

            const size_t n_params = params.size();

            // TODO: Come up with a better way to ensure consistency
            // Perhaps store list/set of free params on setup?
            if(!((n_params + 1) == grads.size())) {
                throw std::runtime_error(
                    "n_params+1=" + std::to_string(n_params+1) + " != grads.size()="
                    + std::to_string(grads.size()) + "; did you set params free/fixed since setup?"
                );
            }

            double loglike_jac = evaluator.loglike_pixel();
            const auto & grad = grads[0];
            size_t n_nonzero = 0;
            for(unsigned int i = 0; i < n_cols; ++i)
            {
                for(unsigned int j = 0; j < n_rows; ++j)
                {
                    n_nonzero += grad.get_value(j, i) != 0;
                }
            }
            if(n_nonzero > 0) {
                errors.push_back("n_nonzero grads[0] entries=" + std::to_string(n_nonzero));
            }

            auto result = _make_evaluator(idx, EvaluatorMode::loglike_image);
            double loglike_img = result.first->loglike_pixel();
            if(loglike_jac != loglike_img) errors.push_back("loglike_jac=" + std::to_string(loglike_jac)
                + "loglike_img=" + std::to_string(loglike_img) + "for exp[" + std::to_string(idx) + "]:");

            const auto & image_current = *(result.second.first);

            for(size_t idx_param = 0; idx_param < n_params; idx_param++) {
                const auto & grad = grads[idx_param + 1];

                auto & param = params[idx_param].get();
                const double value = param.get_value_transformed();
                double diff = value*findiff_frac;
                if(std::abs(diff) < findiff_add) diff = findiff_add;
                double value_new = value + diff;

                auto result_diff = _make_evaluator(idx, EvaluatorMode::image);

                try {
                    param.set_value_transformed(value_new);
                    // If the param is near an upper limit, this might fail: try -diff then
                    if(param.get_value_transformed() != value_new) {
                        diff = -diff;
                        value_new = value + diff;
                        param.set_value_transformed(value_new);
                    }
                } catch(const std::runtime_error & err) {
                    diff = -diff;
                    value_new = value + diff;
                    param.set_value_transformed(value_new);
                    if(param.get_value_transformed() != value_new) {
                        throw std::runtime_error("Couldn't set param=" + param.str()
                            + " to new value_transformed=" + std::to_string(value) + " +/- " 
                            + std::to_string(diff) + "; are limits too restrictive?"
                        );
                    }
                }

                result_diff.first->loglike_pixel();

                size_t n_failed = 0;
                std::vector<double> ratios = {};
                for(unsigned int i = 0; i < n_cols; ++i)
                {
                    for(unsigned int j = 0; j < n_rows; ++j)
                    {
                        double delta = sigma_inv.get_value(j,i)/diff*(
                            result_diff.second.first->get_value(j, i) - image_current.get_value(j, i));
                        double jac = grad.get_value(j, i);

                        auto close = isclose(jac, delta, rtol, atol);
                        if(!close.isclose) {
                            n_failed++;
                            ratios.push_back(jac/delta);
                        }
                    }
                }

                param.set_value_transformed(value);
                if(param.get_value_transformed() != value) {
                    throw std::logic_error("Could not return param=" + param.str() + " to original value="
                        + std::to_string(value));
                }
                if(n_failed > 0) {
                    std::sort(ratios.begin(), ratios.end());
                    double median = ratios[ratios.size()/2];
                    errors.push_back(
                        "Param[" + std::to_string(idx_param) + "]=" + param.str() + " failed for exp["
                        + std::to_string(idx) + "]: n_failed=" + std::to_string(n_failed)
                        + "; median computed/expected=" + std::to_string(median)
                    );
                }
            }
        }
        return errors;
    }

    /**
     * Construct a Model instance from Data, PsfModels and Sources.
     *
     * @param data The data to model.
     * @param psfmodels A vector of PSF models, ordered to match each Observation in data.
     * @param sources A vector of Source models.
     */
    Model(
        std::shared_ptr<const ModelData> data,
        PsfModels & psfmodels,
        Sources & sources
    ) : _data(std::move(data)), _size(_data == nullptr ? 0 : _data->size()) {
        if(_data == nullptr) throw std::invalid_argument("Model data can't be null");
        size_t size = this->size();
        if(psfmodels.size() != size) {
            throw std::invalid_argument("Model psfmodels.size()=" + std::to_string(psfmodels.size())
                + "!= data.size()=" + std::to_string(size));
        }

        _outputs.resize(size);
        _psfmodels.reserve(size);
        size_t i = 0;
        for(auto & psfmodel : psfmodels) {
            if(psfmodel == nullptr) throw std::invalid_argument("Model PsfModels[" + std::to_string(i)
                + "] can't be null");
            _psfcomps.push_back(psfmodel->get_gaussians());
            _psfmodels.push_back(std::move(psfmodel));
            i++;
        }

        _sources.reserve(sources.size());
        i = 0;
        for(auto & source : sources) {
            if(source == nullptr) throw std::invalid_argument("Model Sources[" + std::to_string(i)
                + "] can't be null");
            _sources.push_back(std::move(source));
            i++;
        }

        _gaussians_srcs = this->_get_gaussians_srcs();
        _factors_extra_in.resize(_size);
        _factors_grad_in.resize(_size);
        _map_extra_in.resize(_size);
        _map_grad_in.resize(_size);
    }
};

} // namespace gauss2d::fit

#endif
