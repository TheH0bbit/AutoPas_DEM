/**
 * @file TrapezGrid.h
 * @author R. Penz
 * @date 11/09/21
 */
#pragma once

#include <functional>
#include <numeric>

#include "Objects.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/TrapezGenerator.h"

/**
 * Class describing a Trapez particle grid object.
 */
class TrapezGrid : public Object {
 public:
  /**
   * Constructor.
   * @param velocity
   * @param typeId
   * @param epsilon
   * @param sigma
   * @param mass
   * @param dimensions
   * @param particleSpacing
   * @param bottomLeftCorner
   * @param radius
   * @param young
   * @param poisson
   */
  TrapezGrid(const std::array<double, 3> &velocity, unsigned long typeId, double epsilon, double sigma, double mass,
           const std::array<size_t, 3> &dimensions, double particleSpacing, 
           const std::array<double, 3> &bottomLeftCorner, double radius, double young, double poisson, bool inXDirection)
      : Object(velocity, typeId, epsilon, sigma, mass),
        dimensions(dimensions),
        particleSpacing(particleSpacing),
        bottomLeftCorner(bottomLeftCorner),
        radius(radius),
        young(young),
        poisson(poisson),
        inXDirection(inXDirection) {}

  [[nodiscard]] double getParticleSpacing() const override { return particleSpacing; }

  /**
   * Getter for ParticlesPerDim
   * @return particlePerDim
   */
  [[nodiscard]] const std::array<size_t, 3> &getParticlesPerDim() const { return particlesPerDim; }

  [[nodiscard]] size_t getParticlesTotal() const override {
    if(inXDirection){
      //effective length (base to top) of the trapez 
      size_t effLength = (unsigned long) sqrt(dimensions[1]*dimensions[1]+dimensions[2]*dimensions[2]);
      return dimensions[0] * effLength + effLength * effLength;
    }
    else
    {
      size_t effLength = (unsigned long) sqrt(dimensions[0]*dimensions[0]+dimensions[2]*dimensions[2]);
      return dimensions[1] * effLength + effLength * effLength;
    }
  }

  [[nodiscard]] std::array<double, 3> getBoxMin() const override { 
    if(inXDirection){
      //effective length (base to top) of the trapez 
      size_t effLength = (unsigned long) sqrt(dimensions[1]*dimensions[1]+dimensions[2]*dimensions[2]);
      std::array<double, 3> boxMin = {bottomLeftCorner[0] - effLength, bottomLeftCorner[1] - effLength, bottomLeftCorner[2] - effLength};
      return boxMin;
    }
    else
    {
      size_t effLength = (unsigned long) sqrt(dimensions[0]*dimensions[0]+dimensions[2]*dimensions[2]);
      std::array<double, 3> boxMin = {bottomLeftCorner[0] - effLength, bottomLeftCorner[1] - effLength, bottomLeftCorner[2] - effLength};
      return boxMin;
    }

  }

  [[nodiscard]] std::array<double, 3> getBoxMax() const override {
    if(inXDirection){
      //effective length (base to top) of the trapez 
      size_t effLength = (unsigned long) sqrt(dimensions[1]*dimensions[1]+dimensions[2]*dimensions[2]);
      std::array<double, 3> boxMax = {bottomLeftCorner[0] + dimensions[0] + effLength, bottomLeftCorner[1] + effLength, bottomLeftCorner[2] + effLength};
      return boxMax;
    }
    else
    {
      size_t effLength = (unsigned long) sqrt(dimensions[0]*dimensions[0]+dimensions[2]*dimensions[2]);
      std::array<double, 3> boxMax = {bottomLeftCorner[0] + effLength, bottomLeftCorner[1] + dimensions[1] + effLength, bottomLeftCorner[2] + effLength};
      return boxMax;
    }
  }

  [[nodiscard]] std::string to_string() const override {
    std::ostringstream output;
    output << "to be implemented" << std::endl;
    /*
    std::setw(_valueOffset) << std::left << "particles-per-dimension"
           << ":  " << autopas::utils::ArrayUtils::to_string(particlesPerDim) << std::endl;
    output << std::setw(_valueOffset) << std::left << "particle-spacing"
           << ":  " << particleSpacing << std::endl;
    output << std::setw(_valueOffset) << std::left << "bottomLeftCorner"
           << ":  " << autopas::utils::ArrayUtils::to_string(bottomLeftCorner) << std::endl;
    output << Object::to_string();
    */
    return output.str();
  }

  void generate(autopas::AutoPas<ParticleType> &autopas) const override {
    ParticleType dummyParticle = getDummyParticle(autopas);
    dummyParticle.setRad(radius);
    dummyParticle.setYoung(young);
    dummyParticle.setPoisson(poisson);
    autopasTools::generators::TrapezGenerator::fillWithParticles(
        autopas, dimensions, dummyParticle, inXDirection, {particleSpacing, particleSpacing, particleSpacing}, bottomLeftCorner);
  }

 private:
  std::array<size_t, 3> particlesPerDim = {0,0,0};
  std::array<size_t, 3> dimensions;
  double particleSpacing;
  std::array<double, 3> bottomLeftCorner;
  double radius;
  double young;
  double poisson;
  bool inXDirection;
};
