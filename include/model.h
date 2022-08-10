#ifndef GAUSS2D_FIT_MODEL_H
#define GAUSS2D_FIT_MODEL_H

#include <map>
#include <memory>

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
    typedef Data<T, Image, Mask> ModelData;
    typedef typename ModelData::Observation Observation;
    typedef std::vector<std::shared_ptr<PsfModel>> PsfModels;
    typedef std::vector<std::shared_ptr<Source>> Sources;

private:
    std::shared_ptr<const ModelData> _data;
    std::vector<std::unique_ptr<Evaluator>> _evaluators = {};
    std::vector<std::shared_ptr<Image>> _outputs = {};
    bool _is_setup = false;
    PsfModels _psfmodels;
    Sources _sources;

public:
    std::vector<double> evaluate() {
        if(!_is_setup) throw std::runtime_error("Can't call evaluate before setup_evaluators");
        
        std::vector<double> result(_evaluators.size());
        size_t i = 0;

        for(auto & evaluator : _evaluators) {
            result[i++] = evaluator->loglike_pixel();
        }

        return result;
    }

    std::shared_ptr<const ModelData> get_data() const { return _data; }

    std::unique_ptr<const gauss2d::Gaussians> get_gaussians(const Channel & channel) const override {
        std::vector<std::optional<const gauss2d::Gaussians::Data>> in;
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
        return params;
    }

    PsfModels get_psfmodels() const {
        return _psfmodels;
    }

    Sources get_sources() const {
        return _sources;
    }

    void setup_evaluators(bool save_gradients=false, bool print=false) {
        if(!_is_setup) {
            _evaluators.reserve(_data->size());
        } else {
            // TODO: check if no-op instead of always re-setting up (i.e. store previous save_gradients)
            _evaluators.resize(0);
        }
        std::map<std::reference_wrapper<const Channel>, std::unique_ptr<const Gaussians>> gaussians_src = {};
        for(const Channel & channel : _data->get_channels()) {
            gaussians_src[channel] = this->get_gaussians(channel);
        }
        size_t i = 0;
        for(auto obsit = _data->cbegin(); obsit != _data->cend(); ++obsit) {
            const Observation & observation = *obsit;
            const size_t n_r = observation.get_n_rows();
            const size_t n_c = observation.get_n_cols();

            const auto coordsys = observation.get_image().get_coordsys_ptr_const();
            std::shared_ptr<Image> output = save_gradients ? nullptr : 
                std::make_shared<Image>(n_r, n_c, coordsys);
            _outputs[i] = output;

            auto data_eval = save_gradients ? observation.get_image_ptr_const() : nullptr;
            auto sigma_inv = save_gradients ? observation.get_sigma_inverse_ptr_const() : nullptr;

            std::shared_ptr<const Indices> grad_param_map = save_gradients ? 
                std::make_shared<Indices>(n_r, n_c, coordsys) : nullptr;
            std::shared_ptr<const Image> grad_param_factor = save_gradients ? 
                std::make_shared<Image>(n_r, n_c, coordsys) : nullptr;

            std::shared_ptr<const Indices> extra_param_map = save_gradients ? 
                std::make_shared<Indices>(n_r, n_c, coordsys) : nullptr;
            std::shared_ptr<const Image> extra_param_factor = save_gradients ?
                std::make_shared<Image>(n_r, n_c, coordsys) : nullptr;

            ConvolvedGaussians::Data data {};
            const std::unique_ptr<const gauss2d::Gaussians> psfcomps = _psfmodels[i]->get_gaussians();
            // bconst size_t n_psf_comps = psfcomps->size();
            for(auto source : std::as_const(*gaussians_src[observation.get_channel()])) {
                if(save_gradients)
                {
                    // TODO: finish this
                    // get source param info
                    throw std::logic_error("save_gradients not yet implemented");
                }
                for(auto psfit = psfcomps->cbegin(); psfit < psfcomps->cend(); psfit++) {
                    const std::shared_ptr<const Gaussian> psfgauss = *psfit;
                    data.emplace_back(std::make_shared<ConvolvedGaussian>(
                        source, *psfit
                    ));
                    if(save_gradients)
                    {
                        // TODO: finish this
                        // get psf param info
                    }
                }
            }
            auto gaussians = std::make_shared<const ConvolvedGaussians>(data);

            std::shared_ptr<ImageArray<double, Image>> grads;
            if(save_gradients) {
                typename ImageArray<double, Image>::Data arrays {};
                for(size_t n_cg = 0; n_cg < gaussians->size(); n_cg++) {
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

            _evaluators.emplace_back(
                std::make_unique<Evaluator>(
                    gaussians,
                    coordsys,
                    data_eval,
                    sigma_inv,
                    output,
                    nullptr, // residual,
                    grads
                    /*
                    grad_param_map,
                    grad_param_factor,
                    extra_param_map,
                    extra_param_factor,
                    nullptr
                    */
                )
            );
            i++;
        }
        _is_setup = true;
    }

    std::string str() const override {
        std::string str = "Model(sources=[";
        for(const auto & s : _sources) str += s->str() + ",";
        return str + "], data=" + _data->str() + ")";
    }

    Model(
        std::shared_ptr<const ModelData> data,
        PsfModels & psfmodels,
        Sources & sources
    ) : _data(std::move(data)) {
        if(_data == nullptr) throw std::invalid_argument("Model data can't be null");
        size_t n_exp = _data->size();
        if(psfmodels.size() != n_exp) {
            throw std::invalid_argument("Model psfmodels.size()=" + std::to_string(psfmodels.size())
                + "!= data.size()=" + std::to_string(n_exp));
        }

        _outputs.resize(n_exp);
        _psfmodels.reserve(n_exp);
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