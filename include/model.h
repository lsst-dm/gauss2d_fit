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

namespace gauss2d
{
namespace fit
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

    enum class EvaluatorMode {
        image, loglike, loglike_image, loglike_grad, jacobian
    };

private:
    std::shared_ptr<const ModelData> _data;
    std::vector<std::unique_ptr<Evaluator>> _evaluators = {};
    std::vector<std::shared_ptr<ImageArray<double, Image>>> _grads;
    bool _is_setup = false;
    EvaluatorMode _mode;
    std::vector<std::shared_ptr<Image>> _outputs = {};
    PsfModels _psfmodels;
    size_t _size;
    Sources _sources;

    void _check_obs_idx(size_t idx) const {
        const size_t n_obs = _data->size();
        if(!(idx < n_obs)) throw std::invalid_argument("idx=" + std::to_string(idx) + " !< data->size()="
            + std::to_string(n_obs));
    }

    GaussiansMap _get_gaussians_srcs() const {
        // A map of gaussians for each source by channel
        // The number of params will be the same but integrals will differ
        GaussiansMap gaussians_srcs{};
        size_t n_gaussians_src = 0;

        for(const Channel & channel : _data->get_channels()) {
            size_t n_gaussians_c = 0;
            gaussians_srcs[channel] = std::vector<std::unique_ptr<const Gaussians>>();
            auto & gaussians_src = gaussians_srcs[channel];
            for(const auto & src : get_sources()) {
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

    // Setup an evaluator with the necessary outputs
    std::pair<
        std::unique_ptr<Evaluator>,
        std::pair<std::shared_ptr<Image>, std::shared_ptr<ImageArray<double, Image>>>
    > _make_evaluator(
        const size_t idx_obs,
        const GaussiansMap & gaussians_srcs,
        EvaluatorMode mode=EvaluatorMode::image,
        std::vector<std::shared_ptr<Image>> outputs={},
        std::shared_ptr<Image> residual=nullptr,
        bool print=false
    ) const {
        _check_obs_idx(idx_obs);
        // TODO: implement
        if(mode == EvaluatorMode::loglike_grad) throw std::runtime_error(
            "loglike_grad mode not implemented yet");
        const Observation & observation = _data->at(idx_obs).get();
        const PsfModel & psfmodel = *(_psfmodels[idx_obs]);
        const Channel & channel = observation.get_channel();
        const size_t n_c = observation.get_n_cols();
        const size_t n_r = observation.get_n_rows();

        const size_t n_outputs = outputs.size();
        const bool has_outputs = n_outputs > 0;
        if(has_outputs) {
            if(mode != EvaluatorMode::jacobian) {
                throw std::invalid_argument("outputs not implemented for non-jacobian modes");
            }
        }

        const bool is_mode_image = mode == EvaluatorMode::image;
        const bool is_mode_jacobian = mode == EvaluatorMode::jacobian;
        const bool do_output = is_mode_image || (mode == EvaluatorMode::loglike_image);

        const auto coordsys = observation.get_image().get_coordsys_ptr_const();
        std::shared_ptr<Image> output = !do_output ? nullptr : 
            std::make_shared<Image>(n_r, n_c, coordsys);

        auto data_eval = !is_mode_image ? observation.get_image_ptr_const() : nullptr;
        auto sigma_inv = !is_mode_image ? observation.get_sigma_inverse_ptr_const() : nullptr;

        const auto & gaussians_src = gaussians_srcs.at(channel);
        const std::unique_ptr<const gauss2d::Gaussians> psfcomps = psfmodel.get_gaussians();
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
                //std::as_const(*
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

        std::shared_ptr<extra_param_map> map_extra;
        std::shared_ptr<grad_param_map> map_grad;
        std::shared_ptr<extra_param_factors> factors_extra;
        std::shared_ptr<grad_param_factors> factors_grad;

        if(is_mode_jacobian) {
            auto map_extra_mut = std::make_shared<Indices>(n_gaussians_conv, 2, coordsys);
            auto map_grad_mut = std::make_shared<Indices>
            (n_gaussians_conv, gauss2d::N_PARAMS_GAUSS2D, coordsys);
            auto factors_extra_mut = std::make_shared<Image>(n_gaussians_conv, 3, coordsys);
            auto factors_grad_mut = std::make_shared<Image>(
                n_gaussians_conv, gauss2d::N_PARAMS_GAUSS2D, coordsys);

            map_extra = std::make_shared<extra_param_map>();
            map_grad = std::make_shared<grad_param_map>();

            factors_extra = std::make_shared<extra_param_factors>();
            factors_grad = std::make_shared<grad_param_factors>();

            map_extra->reserve(n_gaussians_conv);
            map_grad->reserve(n_gaussians_conv);
            factors_extra->reserve(n_gaussians_conv);
            factors_grad->reserve(n_gaussians_conv);

            ParameterMap offsets {};
            for(size_t i_psf = 0; i_psf < n_gaussians_psf; ++i_psf)
            {
                for(const auto & src : this->get_sources())
                {
                    src->add_grad_param_map(channel, *map_grad, offsets);
                    src->add_extra_param_map(channel, *map_extra, *map_grad, offsets);

                    src->add_grad_param_factors(channel, *factors_grad);
                    src->add_extra_param_factors(channel, *factors_extra);
                }
            }
            if(map_grad->size() != n_gaussians_conv) {
                throw std::logic_error(
                    "map_grad.size()=" + std::to_string(map_grad->size()) + "!= n_gaussians_conv="
                    + std::to_string(n_gaussians_conv)
                );
            }
            for(size_t row = 0; row < n_gaussians_conv; ++row)
            {
                for(size_t col = 0; col < gauss2d::N_EXTRA_MAP; ++col) {
                    map_extra_mut->set_value_unchecked(row, col, (*map_extra)[row][col]);
                    // std::cerr << map_extra_mut->get_value_unchecked(row, col) << ",";
                }
                // std::cerr << std::endl;
                for(size_t col = 0; col < gauss2d::N_EXTRA_FACTOR; ++col) {
                    factors_extra_mut->set_value_unchecked(row, col, (*factors_extra)[row][col]);
                }
                for(size_t col = 0; col < gauss2d::N_PARAMS_GAUSS2D; ++col) {
                    map_grad_mut->set_value_unchecked(row, col, (*map_grad)[row][col]);
                    factors_grad_mut->set_value_unchecked(row, col, (*factors_grad)[row][col]);
                    // std::cerr << map_grad_mut->get_value_unchecked(row, col) << ",";
                }
                // std::cerr << std::endl;
            }
            map_extra_in = std::move(map_extra_mut);
            map_grad_in = std::move(map_grad_mut);
            factors_extra_in = std::move(factors_extra_mut);
            factors_grad_in = std::move(factors_grad_mut);
        }
        
        std::shared_ptr<ImageArray<double, Image>> grads;
        if(is_mode_jacobian) {
            typename ImageArray<double, Image>::Data arrays {};
            ParamCRefs params_free {};
            auto filter_free = g2f::ParamFilter{false, true, true, true, channel};
            get_parameters_observation_const(params_free, idx_obs, &filter_free);
            params_free = nonconsecutive_unique<ParamBaseCRef>(params_free);
            //for(const auto & param : free_unique) std::cerr << param.get().str() << endl;
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

public:
    void add_extra_param_map(const Channel & channel, extra_param_map & map_extra, const grad_param_map & map_grad, ParameterMap & offsets
        ) const override
    {
        for(auto & source : _sources) source->add_extra_param_map(channel, map_extra, map_grad, offsets);
    }
    void add_extra_param_factors(const Channel & channel, extra_param_factors & factors) const override
    {
        for(auto & source : _sources) source->add_extra_param_factors(channel, factors);
    }
    void add_grad_param_map(const Channel & channel, grad_param_map & map, ParameterMap & offsets
        ) const override
    {
        for(auto & source : _sources) source->add_grad_param_map(channel, map, offsets);
    }
    void add_grad_param_factors(const Channel & channel, grad_param_factors & factors) const override
    {
        for(auto & source : _sources) source->add_grad_param_factors(channel, factors);
    }
    
    std::vector<double> evaluate() {
        if(!_is_setup) throw std::runtime_error("Can't call evaluate before setup_evaluators");

        std::vector<double> result(_evaluators.size());
        size_t i = 0;

        for(auto & evaluator : _evaluators) {
            result[i++] = evaluator->loglike_pixel();
        }

        return result;
    }

    double evaluate(size_t idx) {
        _check_obs_idx(idx);
        return _evaluators[idx]->evaluate();
    }

    std::shared_ptr<const ModelData> get_data() const { return _data; }

    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override {
        std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
        // TODO: Rethink this; it's not sufficient unless sources are single component gaussians
        in.reserve(_sources.size());

        for(auto & source : _sources) {
            in.push_back(source->get_gaussians(channel)->get_data());
        }

        return std::make_unique<gauss2d::Gaussians>(in);
    }

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

    ParamRefs & get_parameters_observation(ParamRefs & params, size_t idx, ParamFilter * filter = nullptr)
    const {
        _check_obs_idx(idx);
        for(auto & source : _sources) source->get_parameters(params, filter);
        _psfmodels[idx]->get_parameters(params, filter);
        _data->at(idx).get().get_parameters(params, filter);
        return params;
    }

    ParamCRefs & get_parameters_observation_const(ParamCRefs & params, size_t idx,
        ParamFilter * filter = nullptr) const {
        _check_obs_idx(idx);
        for(auto & source : _sources) source->get_parameters_const(params, filter);
        _psfmodels[idx]->get_parameters_const(params, filter);
        _data->at(idx).get().get_parameters_const(params, filter);
        return params;
    }

    PsfModels get_psfmodels() const {
        return _psfmodels;
    }

    Sources get_sources() const {
        return _sources;
    }

    // Setup evaluators for all observations for repeated evaluations 
    // (e.g. for fitting with modified params)
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
            
            auto gaussians_srcs = _get_gaussians_srcs();
            size_t i = 0;

            // Apparently this can't go directly in a ternary statement, so here it is
            std::vector<std::shared_ptr<Image>> outputs_obs{};

            for(auto obsit = _data->cbegin(); obsit != _data->cend(); ++obsit) {
                auto result = this->_make_evaluator(
                    i, gaussians_srcs, mode,
                    has_outputs ? outputs[i] : outputs_obs,
                    has_residuals ? residuals[i] : nullptr,
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
                i++;
            }
            _is_setup = true;
            _mode = mode;
        }
    }

    size_t size() const { return _size; }

    std::string str() const override {
        std::string str = "Model(sources=[";
        for(const auto & s : _sources) str += s->str() + ",";
        return str + "], data=" + _data->str() + ")";
    }

    std::vector<std::string> verify_jacobian(
        double findiff_frac=1e-4,
        double findiff_add=1e-4,
        double rtol=1e-3,
        double atol=1e-3
    ) {
        if(_mode != EvaluatorMode::jacobian) this->setup_evaluators(EvaluatorMode::jacobian);

        auto gaussians_srcs = _get_gaussians_srcs();

        const size_t n_exp = _data->size();
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

            auto result = _make_evaluator(idx, gaussians_srcs, EvaluatorMode::loglike_image);
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

                auto result_diff = _make_evaluator(idx, gaussians_srcs, EvaluatorMode::image);

                param.set_value_transformed(value_new);
                // If the param is near an upper limit, this might fail: try -diff then
                if(param.get_value_transformed() != value_new) {
                    diff = -diff;
                    value_new = value + diff;
                    param.set_value_transformed(value_new);
                }

                result_diff.first->loglike_pixel();

                size_t n_failed = 0;
                for(unsigned int i = 0; i < n_cols; ++i)
                {
                    for(unsigned int j = 0; j < n_rows; ++j)
                    {
                        double delta = sigma_inv.get_value(j,i)/diff*(
                            result_diff.second.first->get_value(j, i) - image_current.get_value(j, i));
                        double jac = grad.get_value(j, i);

                        auto close = isclose(jac, delta, rtol, atol);
                        n_failed += !close.isclose;
                    }
                }

                param.set_value_transformed(value);
                if(!(param.get_value_transformed() == value)) {
                    throw std::logic_error("Could not return param=" + param.str() + " to original value="
                        + std::to_string(value));
                }
                if(n_failed > 0) {
                    errors.push_back(
                        "Param[" + std::to_string(idx_param) + "]=" + param.str() + " failed for exp["
                        + std::to_string(idx) + "]: n_failed=" + std::to_string(n_failed)
                    );
                }
            }
        }
        return errors;
    }

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
    }
};

} // namespace fit
} // namespace gauss2d

#endif
