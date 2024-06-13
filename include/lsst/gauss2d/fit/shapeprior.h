#ifndef LSST_GAUSS2D_FIT_SHAPEPRIOR_H
#define LSST_GAUSS2D_FIT_SHAPEPRIOR_H

#include <memory>

#include "lsst/gauss2d/object.h"

#include "parametricellipse.h"
#include "parametricgaussian1d.h"
#include "prior.h"
#include "transforms.h"

namespace lsst::gauss2d::fit {

/**
 * Options for a ShapePrior.
 */
class ShapePriorOptions : public Object {
public:
    static inline const double delta_jacobian_default = 1e-5;
    static inline const double size_maj_floor_default = 1e-3;
    static inline const double axrat_floor_default = 1e-3;

    /**
     * Construct a ShapePriorOptions from values.
     *
     * @param delta_jacobian The small value difference to add to parameter values when computing Jacobians
     * @param size_maj_floor The floor value for the major axis size
     * @param axrat_floor The floor value for the axis ratio
     */
    explicit ShapePriorOptions(double delta_jacobian = delta_jacobian_default,
                               double size_maj_floor = size_maj_floor_default,
                               double axrat_floor = axrat_floor_default);

    bool check_delta_jacobian(double delta_jacobian, bool do_throw = false);
    bool check_size_maj_floor(double size_maj_floor, bool do_throw = false);
    bool check_axrat_floor(double axrat_floor, bool do_throw = false);

    double get_delta_jacobian() const;
    double get_size_maj_floor() const;
    double get_axrat_floor() const;

    void set_delta_jacobian(double delta_jacobian);
    void set_size_maj_floor(double size_maj_floor);
    void set_axrat_floor(double axrat_floor);

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    double _delta_jacobian;
    double _size_maj_floor;
    double _axrat_floor;
};

/**
 * A two-part prior on the shape of a parametric ellipse.
 *
 * @note The size and axis ratio priors are separate and optional.
 * @note The size prior applies to the major-axis size.
 */
class ShapePrior : public Prior {
public:
    /**
     * Construct a ShapePrior from a Parameter and mean_size/std. deviation.
     *
     * @param ellipse The ParametricEllipse to compute a prior for.
     * @param mean_size The mean value of the size prior.
     * @param stddev_size The standard deviation of the size prior.
     */
    explicit ShapePrior(std::shared_ptr<const ParametricEllipse> ellipse,
                        std::shared_ptr<ParametricGaussian1D> prior_size = nullptr,
                        std::shared_ptr<ParametricGaussian1D> prior_axrat = nullptr,
                        std::shared_ptr<ShapePriorOptions> options = nullptr);
    ~ShapePrior();

    PriorEvaluation evaluate(bool calc_jacobians = false, bool normalize_loglike = false) const override;

    std::shared_ptr<ParametricGaussian1D> get_prior_size() const;
    std::shared_ptr<ParametricGaussian1D> get_prior_axrat() const;

    std::vector<double> get_loglike_const_terms() const override;

    void set_prior_size(std::shared_ptr<ParametricGaussian1D> prior_size);
    void set_prior_axrat(std::shared_ptr<ParametricGaussian1D> prior_axrat);

    size_t size() const override;

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override;
    std::string str() const override;

private:
    std::shared_ptr<const ParametricEllipse> _ellipse;
    std::shared_ptr<ParametricGaussian1D> _prior_size;
    std::shared_ptr<ParametricGaussian1D> _prior_axrat;
    std::shared_ptr<ShapePriorOptions> _options;
};
}  // namespace lsst::gauss2d::fit

#endif
