#ifndef GAUSS2D_FIT_MODEL_H
#define GAUSS2D_FIT_MODEL_H

#include <map>
#include <memory>
#include <stdexcept>

#include "gauss2d/evaluate.h"

#include "data.h"
#include "parametricmodel.h"
#include "param_filter.h"
#include "psfmodel.h"
#include "source.h"

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

private:
    std::shared_ptr<const ModelData> _data;
    std::vector<std::unique_ptr<Evaluator>> _evaluators = {};
    std::vector<std::shared_ptr<ImageArray<double, Image>>> _grads;
    bool _is_setup = false;
    bool _is_setup_gradients = false;
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

    std::pair<
        std::unique_ptr<Evaluator>,
        std::pair<std::shared_ptr<Image>, std::shared_ptr<ImageArray<double, Image>>>
    > _make_evaluator(
        const Observation & observation,
        const PsfModel & psfmodel, 
        const GaussiansMap & gaussians_srcs,
        bool save_gradients=false,
        bool print=false
    ) const {
        const Channel & channel = observation.get_channel();
        const size_t n_c = observation.get_n_cols();
        const size_t n_r = observation.get_n_rows();

        const auto coordsys = observation.get_image().get_coordsys_ptr_const();
        std::shared_ptr<Image> output = save_gradients ? nullptr : 
            std::make_shared<Image>(n_r, n_c, coordsys);

        auto data_eval = save_gradients ? observation.get_image_ptr_const() : nullptr;
        auto sigma_inv = save_gradients ? observation.get_sigma_inverse_ptr_const() : nullptr;

        const auto & gaussians_src = gaussians_srcs.at(channel);
        const std::unique_ptr<const gauss2d::Gaussians> psfcomps = psfmodel.get_gaussians();
        const size_t n_gaussians_psf = psfcomps->size();
        size_t n_gaussians_src = 0;
        for(const auto & gaussians : gaussians_src) n_gaussians_src += gaussians->size();
        const size_t n_gaussians_conv = n_gaussians_src*n_gaussians_psf;

        std::shared_ptr<const Indices> map_extra_in, map_grad_in;
        std::shared_ptr<const Image> factors_extra_in, factors_grad_in;

        std::shared_ptr<extra_param_map> map_extra;
        std::shared_ptr<grad_param_map> map_grad;
        std::shared_ptr<extra_param_factors> factors_extra;
        std::shared_ptr<grad_param_factors> factors_grad;

        if(save_gradients) {
            auto map_extra_mut = std::make_shared<Indices>(n_gaussians_conv, 2, coordsys);
            auto map_grad_mut = std::make_shared<Indices>(n_gaussians_conv, gauss2d::N_PARAMS, coordsys);
            auto factors_extra_mut = std::make_shared<Image>(n_gaussians_conv, 3, coordsys);
            auto factors_grad_mut = std::make_shared<Image>(n_gaussians_conv, gauss2d::N_PARAMS, coordsys);

            map_extra = std::make_shared<extra_param_map>();
            map_grad = std::make_shared<grad_param_map>();

            factors_extra = std::make_shared<extra_param_factors>();
            factors_grad = std::make_shared<grad_param_factors>();

            map_extra->reserve(n_gaussians_conv);
            map_grad->reserve(n_gaussians_conv);
            factors_extra->reserve(n_gaussians_conv);
            factors_grad->reserve(n_gaussians_conv);

            ParameterMap offsets {};
            for(const auto & src : this->get_sources())
            {
                for(size_t i_psf = 0; i_psf < n_gaussians_psf; ++i_psf)
                {
                    src->add_grad_param_map(channel, *map_grad, offsets);
                    src->add_extra_param_map(channel, *map_extra, offsets);

                    src->add_grad_param_factors(channel, *factors_grad);
                    src->add_extra_param_factors(channel, *factors_extra);

                }
                std::cerr << "offsets.size()=" << offsets.size() << ","  << map_grad->size() << std::endl;
            }
            if(map_grad->size() != n_gaussians_conv) throw std::logic_error("map_grad.size()=" + std::to_string(map_grad->size())
                + "!= n_gaussians_conv=" + std::to_string(n_gaussians_conv));
            for(size_t row = 0; row < n_gaussians_conv; ++row)
            {
                for(size_t col = 0; col < gauss2d::N_EXTRA_MAP; ++col) {
                    map_extra_mut->set_value_unchecked(row, col, (*map_extra)[row][col]);
                }
                for(size_t col = 0; col < gauss2d::N_EXTRA_FACTOR; ++col) {
                    factors_extra_mut->set_value_unchecked(row, col, (*factors_extra)[row][col]);
                }
                for(size_t col = 0; col < gauss2d::N_PARAMS; ++col) {
                    map_grad_mut->set_value_unchecked(row, col, (*map_grad)[row][col]);
                    factors_grad_mut->set_value_unchecked(row, col, (*factors_grad)[row][col]);
                }
            }
            map_extra_in = std::move(map_extra_mut);
            map_grad_in = std::move(map_grad_mut);
            factors_extra_in = std::move(factors_extra_mut);
            factors_grad_in = std::move(factors_grad_mut);
        }
        
        ConvolvedGaussians::Data data {};

        for(const auto & gaussians : gaussians_src) {
            for(const auto & gaussian : *gaussians) {
                //std::as_const(*
                for(auto psfit = psfcomps->cbegin(); psfit < psfcomps->cend(); psfit++) {
                    data.emplace_back(std::make_shared<ConvolvedGaussian>(gaussian, *psfit));
                }
            }
        }
        auto gaussians = std::make_shared<const ConvolvedGaussians>(data);
        if(gaussians->size() != n_gaussians_conv) {
            throw std::logic_error("gaussians.size()=" + std::to_string(gaussians->size())
                + "!= n_gaussians_conv=" + std::to_string(n_gaussians_conv) + "= n_gaussians_src="
                + std::to_string(n_gaussians_src) + "= n_gaussians_psf=" + std::to_string(n_gaussians_psf)
                + "(save_gradients=" + std::to_string(save_gradients)
            );
        }

        std::shared_ptr<ImageArray<double, Image>> grads;
        if(save_gradients) {
            typename ImageArray<double, Image>::Data arrays {};
            ParamCRefs free {};
            auto filter_free = g2f::ParamFilter{false, true};
            this->get_parameters_const(free, &filter_free);
            const size_t n_free = free.size() + 1;

            for(size_t idx_img = 0; idx_img < n_free; ++idx_img) {
                arrays.emplace_back(
                    std::make_shared<Image>(n_r, n_c, coordsys)
                );
            }
            grads = std::make_shared<ImageArray<double, Image>>(&arrays);
        } else {
            grads = nullptr;
        }
        if(print) {
            std::cout << observation.str() << std::endl;
            std::cout << gaussians->str() << std::endl;
        }

        // Prevent nulling of grads/output by moving in make_unique<Evaluator>
        auto grads2 = grads;
        auto output2 = output;

        return {
            std::make_unique<Evaluator>(
                gaussians,
                coordsys,
                data_eval,
                sigma_inv,
                output,
                nullptr, // residual,
                grads,
                map_grad_in,
                factors_grad_in,
                map_extra_in,
                factors_extra_in
                //, background
            ),
            {std::move(output2), std::move(grads2)}
        };
    }

public:
    void add_extra_param_map(const Channel & channel, extra_param_map & map, ParameterMap & offsets
        ) const override
    {
        for(auto & source : _sources) source->add_extra_param_map(channel, map, offsets);
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
        _data[idx]->get_parameters(params, filter);
        return params;
    }

    ParamCRefs & get_parameters_observation_const(ParamCRefs & params, size_t idx,
        ParamFilter * filter = nullptr) const {
        _check_obs_idx(idx);
        for(auto & source : _sources) source->get_parameters_const(params, filter);
        _psfmodels[idx]->get_parameters_const(params, filter);
        _data[idx]->get_parameters(params, filter);
        return params;
    }

    PsfModels get_psfmodels() const {
        return _psfmodels;
    }

    Sources get_sources() const {
        return _sources;
    }

    void setup_evaluators(bool save_gradients=false, bool print=false) {
        if(!_is_setup || (!save_gradients == _is_setup_gradients)) {
            _evaluators.resize(0);
            _grads.resize(0);
            _outputs.resize(0);
            _evaluators.reserve(size());
            if(save_gradients) {
                _grads.reserve(size());
            } else {
                _outputs.reserve(size());
            }
            
            const auto & channels = _data->get_channels();
            
            auto gaussians_srcs = _get_gaussians_srcs();
            size_t i = 0;
            for(auto obsit = _data->cbegin(); obsit != _data->cend(); ++obsit) {
                auto result = this->_make_evaluator(*obsit, *_psfmodels[i], gaussians_srcs, save_gradients);
                _evaluators.emplace_back(std::move(result.first));
                if(save_gradients) {
                    _grads.emplace_back(std::move(result.second.second));
                } else {
                    _outputs.emplace_back(std::move(result.second.first));
                }
                i++;
            }
            _is_setup = true;
            _is_setup_gradients = save_gradients;
        }
    }

    size_t size() const { return _size; }

    std::string str() const override {
        std::string str = "Model(sources=[";
        for(const auto & s : _sources) str += s->str() + ",";
        return str + "], data=" + _data->str() + ")";
    }

    bool verify_gradients() {
        if(!this->_is_setup_gradients) this->setup_evaluators(true);

        auto gaussians_srcs = _get_gaussians_srcs();

        const size_t n_exp = _data->size();
        ParamFilter filter{false, true};
        for(size_t idx; idx < n_exp; ++idx)
        {
            const auto & evaluator = _evaluators->at(idx);
            const auto & grad = _grads.at(idx);
            const auto & observation = _data->at(idx);

            const Channel & channel = observation.get_channel();
            const size_t n_c = observation.get_n_cols();
            const size_t n_r = observation.get_n_rows();

            filter.channel = channel;

            ParamRefs params{};
            this->_get_parameters_observation(params, idx, &filter);

            const size_t n_params = params.size();

            // TODO: Come up with a better way to ensure consistency
            // Perhaps store list/set of free params on setup?
            if(!(n_params == _grads.size())) {
                throw std::runtime_error(
                    "n_params=" + std::to_string(n_params) + " != _grads.size()=" + _grads.size()
                    + "; did you set params free/fixed since setup?");
            }

            evaluator->evaluate();

            auto result = _make_evaluator(observation, *(_psfmodels[idx]), gaussians_srcs);
            result.first->evaluate();

            const auto & image_current = result.second;

            for(size_t idx_param; idx_param < n_params; idx_param++) {
                const auto & param = params[idx_param].get();
                const double value = param.get_value_transformed();
            }
        }
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