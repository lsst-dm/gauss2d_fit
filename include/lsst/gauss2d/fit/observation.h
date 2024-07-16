#ifndef LSST_GAUSS2D_FIT_OBSERVATION_H
#define LSST_GAUSS2D_FIT_OBSERVATION_H

#include <memory>
#include <stdexcept>

#include "lsst/gauss2d/image.h"
#include "lsst/gauss2d/type_name.h"

#include "channel.h"
#include "param_defs.h"
#include "parametric.h"
#include "util.h"

namespace lsst::gauss2d::fit {

/**
 * @brief An observed single-channel image with an associated variance and mask.
 *
 * An Observation is an Image in a single Channel with a corresponding variance
 * and mask. The variance is provided as inverse sigma to facilitate Model
 * evaluation. Observation Image instances may represent a single physical
 * exposure or more generic data, such as processed data like coadded/stacked
 * images.
 *
 * @note The Mask convention is that of an "inverse" mask. Pixels with
 *       positive Mask valuese used, while values <= 0 are ignored.
 *
 * @tparam T The type of the Image (usually float or double).
 * @tparam Image The class of image used in Data.
 * @tparam Indices The class of index map used in Data.
 * @tparam Mask The class of mask f.
 */
template <typename T, typename I, typename M>
class Observation : public Parametric {
public:
    using Image = lsst::gauss2d::Image<T, I>;
    using Mask = lsst::gauss2d::Image<bool, M>;

    /**
     * Construct an Observation instance.
     *
     * @param image The Image to assign to _image.
     * @param sigma_inv The Image to assign to _sigma_inv. Must have identical dimensions as image.
     * @param mask_inv The mask Image to assign to _mask.
     * @param channel The channel of every Observation.
     */
    explicit Observation(std::shared_ptr<Image> image, std::shared_ptr<Image> sigma_inv,
                         std::shared_ptr<Mask> mask_inv, const Channel &channel = Channel::NONE())
            : _image(std::move(image)),
              _sigma_inv(std::move(sigma_inv)),
              _mask_inv(std::move(mask_inv)),
              _channel(channel) {
        if ((_image == nullptr) || (_sigma_inv == nullptr) || (_mask_inv == nullptr)) {
            throw std::invalid_argument("Must supply non-null image, variance and mask");
        }
        std::string msg;
        bool passed = images_compatible<T, I, T, I>(*_image, *_sigma_inv, true, &msg);
        passed &= images_compatible<T, I, bool, M>(*_image, *_mask_inv, true, &msg);
        if (passed != (msg.empty()))
            throw std::logic_error("Observation images_compatible=" + std::to_string(passed)
                                   + " != msg != '' (=" + msg + ")");
        if (!passed) throw std::invalid_argument("image/variance/mask incompatible: " + msg);
    }

    /// Get this->_channel
    const Channel &get_channel() const { return _channel; }
    /// Get this->_image
    std::shared_ptr<const Image> get_image_ptr_const() const { return _image; }
    /// Get this->_sigma_inv
    std::shared_ptr<const Image> get_sigma_inverse_ptr_const() const { return _sigma_inv; }

    /// Get the number of columns in this->_image
    size_t get_n_cols() const { return _image->get_n_cols(); }
    /// Get the number of rows in this->_image
    size_t get_n_rows() const { return _image->get_n_rows(); }

    /// Get a ref to this->_image
    Image &get_image() const { return *_image; }
    /// Get a ref to this->_mask
    Mask &get_mask_inverse() const { return *_mask_inv; }
    /// Get a ref to this->sigma_inv
    Image &get_sigma_inverse() const { return *_sigma_inv; }

    ParamRefs &get_parameters(ParamRefs &params, ParamFilter *filter = nullptr) const override {
        return params;
    }
    ParamCRefs &get_parameters_const(ParamCRefs &params, ParamFilter *filter = nullptr) const override {
        return params;
    }

    std::string repr(bool name_keywords = false,
                     std::string_view namespace_separator = Object::CC_NAMESPACE_SEPARATOR) const override {
        return type_name_str<Observation>(false, namespace_separator) + ")" + (name_keywords ? "image=" : "")
               + repr_ptr(_image.get(), name_keywords, namespace_separator) + ", "
               + (name_keywords ? "sigma_inv=" : "")
               + repr_ptr(_sigma_inv.get(), name_keywords, namespace_separator) + ", "
               + (name_keywords ? "mask_inv=" : "")
               + repr_ptr(_mask_inv.get(), name_keywords, namespace_separator) + ", "
               + (name_keywords ? "channel=" : "") + _channel.repr(name_keywords, namespace_separator) + ")";
    }
    std::string str() const override {
        return type_name_str<Observation>(true) + "(image=" + str_ptr(_image.get())
               + ", sigma_inv=" + str_ptr(_sigma_inv.get()) + ", mask_inv=" + str_ptr(_mask_inv.get())
               + ", channel=" + _channel.str() + ")";
    }

    bool operator==(const Observation &other) const {
        return ((this->get_image() == other.get_image())
                && (this->get_mask_inverse() == other.get_mask_inverse())
                && (this->get_sigma_inverse() == other.get_sigma_inverse()));
    }

private:
    std::shared_ptr<Image> _image;
    std::shared_ptr<Image> _sigma_inv;
    std::shared_ptr<Mask> _mask_inv;
    const Channel &_channel;
};

}  // namespace lsst::gauss2d::fit

#endif
