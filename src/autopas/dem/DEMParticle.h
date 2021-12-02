/**
 * @file DEMParticle.h
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
class DEMParticle final : public Particle {
 public:
  DEMParticle() = default;

  /**
   * Constructor of dem particle with initialization of additional parameters.
   * @param pos Position of the particle.
   * @param v Velocitiy of the particle.
   * @param particleId Id of the particle.
   * @param rad Radius of the particle.
   * @param poisson Poisson ratio of the particle.
   * @param young Young modulus of the particle.
   * 
   */
  explicit DEMParticle(std::array<floatType, 3> pos, std::array<floatType, 3> v, unsigned long particleId,
                      floatType rad = 0.0, floatType poisson = 0.0, floatType young = 0.0)
      : Particle(pos, v, particleId), _rad(rad), _poisson(poisson), _young(young) {}

  ~DEMParticle() final = default;

  /**
   * Enums used as ids for accessing and creating a dynamically sized SoA.
   */
  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, rad, poisson, young, mass, typeId, ownershipState };

  /**
   * The type for the SoA storage.
   *
   */
  using SoAArraysType =
      typename autopas::utils::SoAType<DEMParticle<floatType> *, size_t /*id*/, floatType /*x*/, floatType /*y*/,
                                       floatType /*z*/, floatType /*fx*/, floatType /*fy*/, floatType /*fz*/,
                                       floatType /*rad*/, floatType /*poisson*/, floatType /*young*/, floatType /*mass*/,
                                       size_t /*typeId*/, OwnershipState /*ownershipState*/>::Type;

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
    } else if constexpr (attribute == AttributeNames::mass) {
      return getMass();
    } else if constexpr (attribute == AttributeNames::typeId) {
      return getTypeId();
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      utils::ExceptionHandler::exception("DEMParticle::get() unknown attribute {}", attribute);
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
      setRad(value);
    } else if constexpr (attribute == AttributeNames::poisson) {
      setPoisson(value);
    } else if constexpr (attribute == AttributeNames::young) {
      setYoung(value);
    } else if constexpr (attribute == AttributeNames::mass) {
      setMass(value);
    } else if constexpr (attribute == AttributeNames::typeId) {
      setTypeId(value);
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      utils::ExceptionHandler::exception("DEMParticle::set() unknown attribute {}", attribute);
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

  /**
   * Get the Mass of a Particle.
   * @return
   */
  [[nodiscard]] floatType getMass() const { return _mass; }

  /**
   * Set the Mass of a Particle.
   * @param young
   */
  void setMass(floatType mass) { _mass = mass; }

  /**
   * Get TypeId.
   * @return
   */
  [[nodiscard]] size_t getTypeId() const { return _typeId; }

  /**
   * Set the type id of the Molecule.
   * @param typeId
   */
  void setTypeId(size_t typeId) { _typeId = typeId; }


 private:
    /**
   * Particle radius.
   */
  floatType _rad = 0.0;

  /**
   * Particle Poisson ratio.
   */
  floatType _poisson = 0.0;  
  
  /**
   * Particle Young modulus.
   */
  floatType _young = 0.0;

  /**
   * Particle Mass.
   */
  floatType _mass = 1.0;

  /**
   * TypeId, needed for compability
   */
  size_t _typeId = 0;

  /**
   * Old Force of the particle experiences as 3D vector.
   */
  std::array<double, 3> _oldF = {0., 0., 0.};
};

}  // namespace autopas
