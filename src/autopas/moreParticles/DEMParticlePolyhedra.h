/**
 * @file DEMParticlePolyhedra.h
 *
 * @date 07.07.2021
 * @author R. Penz
 */

#pragma once

#include <vector>

#include "autopas/particles/Particle.h"

namespace autopas {

/**
 * Particle class for the DEMFunctor.
 */
template <typename floatType = double>
class DEMParticlePolyhedra final : public Particle {
 public:
  DEMParticlePolyhedra() = default;

  /**
   * @param pos Position of the particle.
   * @param v Velocitiy of the particle.
   * @param particleId Id of the particle.
   * @param rad Radius of the particle.
   * @param poisson Poisson ratio of the particle.
   * @param young Young modulus of the particle.
   * 
   * 
   * 
   * 
   */
  explicit DEMParticlePolyhedra(std::array<floatType, 3> pos, std::array<floatType, 3> v, unsigned long particleId,
                      double rad = 0.0, double poisson = 0.0, double young = 0.0)
      : Particle(pos, v, particleId), _rad(rad), _poisson(poisson), _young(young) {}

  ~DEMParticlePolyhedra() final = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, rad, poisson, young, ownershipState };

  /**
   * The type for the SoA storage.
   *
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<MoleculeLJ<floatType> *, size_t /*id*/, floatType /*x*/, floatType /*y*/,
                                       floatType /*z*/, floatType /*fx*/, floatType /*fy*/, floatType /*fz*/,
                                       floatType /*rad*/, floatType /*poisson*/, floatType /*young*/, OwnershipState /*ownershipState*/>::Type;

  /**
   * Getter, which allows access to an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @return Value of the requested attribute.
   * @note The value of owned is return as floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    if constexpr (attribute == AttributeNames::ptr) {
      return this;
    } else if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::rad) {
      return getRad();
    } else if constexpr (attribute == AttributeNames::poisson) {
      return getPoisson();
    } else if constexpr (attribute == AttributeNames::young) {
      return getYoung();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      utils::ExceptionHandler::exception("DEMParticlePolyhedra::get() unknown attribute {}", attribute);
    }
  }

  /**
   * Setter, which allows set an attribute using the corresponding attribute name (defined in AttributeNames).
   * @tparam attribute Attribute name.
   * @param value New value of the requested attribute.
   * @note The value of owned is extracted from a floating point number (true = 1.0, false = 0.0).
   */
  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::id) {
      setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      _r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      _r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      _r[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      _f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      _f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      _f[2] = value;
    } else if constexpr (attribute == AttributeNames::rad) {
      setRad();
    } else if constexpr (attribute == AttributeNames::poisson) {
      setPoisson();
    } else if constexpr (attribute == AttributeNames::young) {
      setYoung();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      utils::ExceptionHandler::exception("DEMParticlePolyhedra::set() unknown attribute {}", attribute);
    }
  }

  /**
   * Get the old force.
   * @return
   */
  [[nodiscard]] std::array<double, 3> getOldf() const { return _oldF; }

  /**
   * Set old force.
   * @param oldForce
   */
  void setOldF(const std::array<double, 3> &oldForce) { _oldF = oldForce; }

  /**
   * Get the Radius of a Particle.
   * @return
   */
  [[nodiscard]] floatType getRad() const { return _rad; }

  /**
   * Set the Radius of a Particle.
   * @param rad
   */
  void setRad(floatType rad) { _rad = rad; }

  /**
   * Get the Poisson ratio of a Particle.
   * @return
   */
  [[nodiscard]] floatType getPoisson() const { return _poisson; }

  /**
   * Set the Poisson ratio of a Particle.
   * @param poisson
   */
  void setPoisson(floatType poisson) { _poisson = poisson; }

  /**
   * Get the Young Modulus of a Particle.
   * @return
   */
  [[nodiscard]] floatType getYoung() const { return _young; }

  /**
   * Set the Young Modulus of a Particle.
   * @param young
   */
  void setYoung(floatType young) { _young = young; }


 private:
    /**
   * Particle radius.
   */
  double _rad = 0.0;

  /**
   * Particle Poisson ratio.
   */
  double _poisson = 0.0;  
  
  /**
   * Particle Young modulus.
   */
  double _young = 0.0;

  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> _oldF = {0., 0., 0.};
};

}  // namespace autopas
