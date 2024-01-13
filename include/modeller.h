#ifndef GAUSS2D_FIT_MODELLER_H
#define GAUSS2D_FIT_MODELLER_H

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "gauss2d/evaluate.h"
#include "gauss2d/image.h"

#include "data.h"
#include "model.h"
#include "param_filter.h"
#include "parameters.h"
#include "prior.h"
#include "psfmodel.h"
#include "source.h"
#include "util.h"

namespace gauss2d::fit {

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
class Modeller : public Parametric {
public:
    typedef Model<T, Image, Indices, Mask> Model;
/*
    typedef GaussianEvaluator<T, Image, Indices> Evaluator;
    typedef std::map<std::reference_wrapper<const Channel>, std::vector<std::unique_ptr<const Gaussians>>>
            GaussiansMap;
    typedef Data<T, Image, Mask> ModelData;
    typedef typename ModelData::Observation Observation;
    typedef std::vector<std::shared_ptr<PsfModel>> PsfModels;
    typedef std::vector<std::shared_ptr<Source>> Sources;
    typedef std::vector<std::shared_ptr<Prior>> Priors;
*/
    /**
     * Construct a Modeller instance.
     */
    Modeller() {
    }

private:

public:
    void fit_model(Model & model, int idx_obs) {

    }

    void fit_model_linear(Model & model, int idx_obs) {

    }
};

}  // namespace gauss2d::fit

#endif  // GAUSS2DFIT_MODELLER_H
