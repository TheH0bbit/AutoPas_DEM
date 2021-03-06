/**
 * @file TimeDiscretization.h
 * @author N. Fottner
 * @date 13/05/19
 */
#pragma once
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {

/**
 * Calculate and update the position for every particle using the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void calculatePositions(AutoPasTemplate &autopas, const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                        const double deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto v = iter->getV();
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
    auto f = iter->getF();
    iter->setOldF(f);
    iter->setF({0., 0., 0.});
    v = mulScalar(v, deltaT);
    f = mulScalar(f, (deltaT * deltaT / (2 * m)));
    auto newR = add(v, f);
    iter->addR(newR);
  }
}


/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void calculateVelocities(AutoPasTemplate &autopas, const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                         const double deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
    auto force = iter->getF();
    auto oldForce = iter->getOldf();
    auto newV = mulScalar((add(force, oldForce)), deltaT / (2 * m));
    iter->addV(newV);
  }
}

/**
 * Calculate and update the position for every DEM particle using the Störmer-Verlet Algorithm.
 * @param autopas
 * @param deltaT time step width
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void calculatePositionsDEM(AutoPasTemplate &autopas, const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                        const double deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto v = iter->getV();
    auto m = iter->getMass();
    auto f = iter->getF();
    iter->setOldF(f);
    iter->setF({0., 0., 0.});
    if(m < 1000){
      v = mulScalar(v, deltaT);
      f = mulScalar(f, (deltaT * deltaT / (2 * m)));
      auto newR = add(v, f);
    iter->addR(newR);
    }
    
  }
}


/**
 * Calculate and update the velocity for every DEM particle using the the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void calculateVelocitiesDEM(AutoPasTemplate &autopas, const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                         const double deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto m = iter->getMass();
    if(m < 1000){
      auto force = iter->getF();
      auto oldForce = iter->getOldf();
      auto newV = mulScalar((add(force, oldForce)), deltaT / (2 * m));
    iter->addV(newV);
    }
  }
}

}  // namespace TimeDiscretization
