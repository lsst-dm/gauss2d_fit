#ifndef GAUSS2D_FIT_OBSERVATION_H
#define GAUSS2D_FIT_OBSERVATION_H

#include <memory>
#include <stdexcept>

#include "gauss2d/image.h"

#include "channel.h"
#include "param_defs.h"
#include "parametric.h"
#include "util.h"

namespace gauss2d
{
namespace fit
{

template <typename T, typename I, typename M>
class Observation : public Parametric
{
public:
    using Image = gauss2d::Image<T, I>;
    using Mask = gauss2d::Image<bool, M>;
private:
    std::shared_ptr<Image> _image;
    std::shared_ptr<Image> _sigma_inv;
    std::shared_ptr<Mask> _mask_inv;
    const Channel & _channel;

public:
    const Channel & get_channel() const { return _channel; }

    std::shared_ptr<const Image> get_image_ptr_const() const { return _image; }
    std::shared_ptr<const Image> get_sigma_inverse_ptr_const() const { return _sigma_inv; }

    size_t get_n_cols() const { return _image->get_n_cols(); }
    size_t get_n_rows() const { return _image->get_n_rows(); }

    Image & get_image() const { return *_image; }
    Mask & get_mask_inverse() const { return *_mask_inv; }
    Image & get_sigma_inverse() const { return *_sigma_inv; }

    ParamRefs & get_parameters(ParamRefs & params, ParamFilter * filter = nullptr) const override {
        return params;
    }
    ParamCRefs & get_parameters_const(ParamCRefs & params, ParamFilter * filter = nullptr) const override {
        return params;
    }

    std::string str() const override {
        return "Observation(image=" + str_ptr(_image.get())
            + "," + str_ptr(_sigma_inv.get())
            + "," + str_ptr(_mask_inv.get())
            + "," + _channel.str()
            + ")";
    }

    const bool operator == (const Observation &other) const {
        return (
            (this->get_image() == other.get_image())
            && (this->get_mask_inverse() == other.get_mask_inverse())
            && (this->get_sigma_inverse() == other.get_sigma_inverse())
        );
    }

    Observation(
        std::shared_ptr<Image> image,
        std::shared_ptr<Image> sigma_inv,
        std::shared_ptr<Mask> mask_inv,
        const Channel & channel = Channel::NONE()
    ) :
        _image(std::move(image)), _sigma_inv(std::move(sigma_inv)), _mask_inv(std::move(mask_inv)),
        _channel(channel)
    {
        if((_image == nullptr) || (_sigma_inv == nullptr) || (_mask_inv == nullptr))
        {
            throw std::invalid_argument("Must supply non-null image, variance and mask");
        }
        std::string msg = "";
        bool passed = images_compatible<T, I, T, I>(*_image, *_sigma_inv, &msg);
        passed &= images_compatible<T, I, bool, M>(*_image, *_mask_inv, &msg);
        if(passed != (msg == "")) throw std::logic_error(
            "Observation images_compatible=" + std::to_string(passed) + " != msg != '' (=" + msg + ")");
        if(!passed) throw std::invalid_argument("image/variance/mask incompatible: " + msg);
    }
};

} // namespace fit
} // namespace gauss2d

#endif