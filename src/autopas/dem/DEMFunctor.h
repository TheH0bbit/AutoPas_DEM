/**
 * @file DEMFunctor.h
 *
 * @date 01.07.21
 * @author R. Penz
 */

#pragma once

#include <array>
#include <cmath>

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * A functor to handle interactions between two DEM particles.
 * @tparam Particle The type of particle.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */

template <class Particle, FunctorN3Modes useNewton3 = FunctorN3Modes::Both, bool relevantForTuning = true>
class DEMFunctor
    : public Functor<Particle,
                     DEMFunctor<Particle, useNewton3, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  DEMFunctor() = delete;

  /**
   * Internal, actual constructor.
   * @param cutoff
   */
  explicit DEMFunctor(double cutoff)

      : Functor<Particle, DEMFunctor<Particle, useNewton3, relevantForTuning>>(cutoff) {}

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final { return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both; }

  bool allowsNonNewton3() final {
    return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
  }
        
/**
     * Hertz elastic Solution
     * F = 4/3*E*sqrt(R)*sqrt(delta^3)
     * R = R1*R2/(R1+R2)
     * 1/E = (1-v1^2)/E1 + (1-v2^2)/E2
     * 
     * 
     * needed Parameters: 
     * E1, E2 - respective Young's modulus, need to specify
     * v1, v2 - respective Poisson ratio, need to specify
     * R - Radii
     * delta - approach distance
     */

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    
    auto dr = utils::ArrayMath::sub(i.getR(), j.getR()); //distance between ParticleCenters
    double penDepth = i.getRad()+j.getRad()-utils::ArrayMath::L2Norm(dr);

    if(penDepth <= 0) return; //return if Particles dont intersect
    double e1 = (1-pow(i.getPoisson(), 2))/i.getYoung();
    double e2 = (1-pow(j.getPoisson(), 2))/j.getYoung();
    double e = 1/(e1+e2);

    double r = i.getRad()*j.getRad()/(i.getRad()+j.getRad());

    //calculate Force and ForceVector
    double f = 4/3*e*sqrt(r)*pow(penDepth, 3./2.);
    auto vecf = utils::ArrayMath::mulScalar (utils::ArrayMath::normalize(dr),f);
    i.addF(vecf);
    if(newton3)
    {
      j.subF(vecf);
    }
  }

  void SoAFunctorCalc(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2,  bool single, bool newton3)
  {
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;
    size_t shift = 0;

    const auto *const __restrict id1ptr = soa1.template begin<Particle::AttributeNames::id>();
    const auto *const __restrict id2ptr = soa2.template begin<Particle::AttributeNames::id>();

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    SoAFloatPrecision *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    SoAFloatPrecision *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict rad1ptr = soa1.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson1ptr = soa1.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young1ptr = soa1.template begin<Particle::AttributeNames::young>();
    const auto *const __restrict rad2ptr = soa2.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson2ptr = soa2.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young2ptr = soa2.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownership1ptr = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownership2ptr = soa2.template begin<Particle::AttributeNames::ownershipState>();
    
    for(unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      if(single){shift++;}  //increase shift by 1 if single View
      const auto ownedStateI = ownership1ptr[i];
      if(ownedStateI == OwnershipState::dummy) {return;}

      //accumulating Force directions
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;
    
// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
// shift for SoAFunctorSingle
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
     for (unsigned int j = shift; j < soa2.getNumParticles(); ++j) {
        const auto ownedStateJ = ownership2ptr[j];

        const SoAFloatPrecision drx = x1ptr[i] - x2ptr[j];
        const SoAFloatPrecision dry = y1ptr[i] - y2ptr[j];
        const SoAFloatPrecision drz = z1ptr[i] - z2ptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const SoAFloatPrecision dr = sqrt(dr2);

        const SoAFloatPrecision radI = rad1ptr[i];
        const SoAFloatPrecision radJ = rad2ptr[j];
        const SoAFloatPrecision rad = radI + radJ;

        // Mask away if particles arent intersecting or if j is dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy;

        const SoAFloatPrecision poissonI = poisson1ptr[i];
        const SoAFloatPrecision youngI = young1ptr[i];
        const SoAFloatPrecision poissonJ = poisson2ptr[j];
        const SoAFloatPrecision youngJ = young2ptr[j];

        const SoAFloatPrecision poissonI2 = poissonI * poissonI;
        const SoAFloatPrecision poissonJ2 = poissonJ * poissonJ;
        const SoAFloatPrecision e1 = (1. - poissonI2)/youngI;
        const SoAFloatPrecision e2 = (1. - poissonJ2)/youngJ;
        const SoAFloatPrecision esum = e1 + e2;
        const SoAFloatPrecision e = 1. / esum;

        const SoAFloatPrecision penDepth = rad - dr;

        const SoAFloatPrecision r = radI * radJ / rad;
        const SoAFloatPrecision fac = mask ? 4/3 * e * sqrt(r) * pow(penDepth, 3./2.) / dr : 0.;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        if(newton3){
          fx2ptr[j] -=fx;
          fy2ptr[j] -=fy;
          fz2ptr[j] -=fz;
        }
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

  /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) final {

    SoAFunctorCalc(soa, soa, true, newton3);

  }

  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2) {
    SoAFunctorCalc(soa1, soa2, false, newton3);
  }

 public:
  // clang-format off
  /**
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   */
  // clang-format on
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
    if (soa.getNumParticles() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 12>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::rad,    Particle::AttributeNames::poisson, 
        Particle::AttributeNames::young,  Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::rad,    Particle::AttributeNames::poisson, 
        Particle::AttributeNames::young,  Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   * Get the number of flops used per kernel call. This should count the
   * floating point operations needed for two particles that lie within a cutoff
   * radius.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply
    // scale) sum Forces: 6 (forces) kernel total = 12 + 6 = 18
    return 18ul;
  }

  void initTraversal() final {
    _postProcessed = false;
  }

  void endTraversal(bool newton3) final {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
  }

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    auto *const __restrict radptr = soa.template begin<Particle::AttributeNames::rad>();
    auto *const __restrict poissonptr = soa.template begin<Particle::AttributeNames::poisson>();
    auto *const __restrict youngptr = soa.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;
    
    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    // checks whether particle i is owned.
    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == OwnershipState::dummy) {
      return;
    }

    // this is a magic number, that should correspond to at least
    // vectorization width*N have testet multiple sizes:
    // 4: does not give a speedup, slower than original AoSFunctor
    // 8: small speedup compared to AoS
    // 12: highest speedup compared to Aos
    // 16: smaller speedup
    // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
    // use a multiple of 8 for avx
    constexpr size_t vecsize = 16;
#else
    // for everything else 12 is faster
    constexpr size_t vecsize = 12;
#endif
    size_t joff = 0;

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, radtmp, poissontmp, youngtmp, xArr, yArr, zArr, radArr, poissonArr, youngArr, fxArr, fyArr, fzArr;
      alignas(64) std::array<OwnershipState, vecsize> ownedStateArr{};
      // broadcast of the position of particle i
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];
        radtmp[tmpj] = radptr[indexFirst];
        poissontmp[tmpj] = poissonptr[indexFirst];
        youngtmp[tmpj] = youngptr[indexFirst];
      }
      // loop over the verlet list from 0 to x*vecsize
      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
        // in each iteration we calculate the interactions of particle i with
        // vecsize particles in the neighborlist of particle i starting at
        // particle joff

        // gather position of particle j
#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xArr[tmpj] = xptr[neighborListPtr[joff + tmpj]];
          yArr[tmpj] = yptr[neighborListPtr[joff + tmpj]];
          zArr[tmpj] = zptr[neighborListPtr[joff + tmpj]];
          radArr[tmpj] = radptr[neighborListPtr[joff + tmpj]];
          poissonArr[tmpj] = poissonptr[neighborListPtr[joff + tmpj]];
          youngArr[tmpj] = youngptr[neighborListPtr[joff + tmpj]];
          ownedStateArr[tmpj] = ownedStatePtr[neighborListPtr[joff + tmpj]];
        }
        // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc) safelen(vecsize)
        for (size_t j = 0; j < vecsize; j++) {

          const auto ownedStateJ = ownedStateArr[j];

          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
          const SoAFloatPrecision dr = sqrt(dr2);

          const SoAFloatPrecision radI = radtmp[j];
          const SoAFloatPrecision radJ = radArr[j];
          const SoAFloatPrecision rad = radI + radJ;

          // Mask away if particles arent intersecting or if j is dummy.
          const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy and dr != 0;

          const SoAFloatPrecision poissonI = poissontmp[j];
          const SoAFloatPrecision youngI = youngtmp[j];
          const SoAFloatPrecision poissonJ = poissonArr[j];
          const SoAFloatPrecision youngJ = youngArr[j];

          const SoAFloatPrecision poissonI2 = poissonI * poissonI;
          const SoAFloatPrecision poissonJ2 = poissonJ * poissonJ;
          const SoAFloatPrecision e1 = (1. - poissonI2)/youngI;
          const SoAFloatPrecision e2 = (1. - poissonJ2)/youngJ;
          const SoAFloatPrecision esum = e1 + e2;
          const SoAFloatPrecision e = 1. / esum;

          const SoAFloatPrecision penDepth = rad - dr;

          const SoAFloatPrecision r = radI * radJ / rad;
          const SoAFloatPrecision fac = mask ? 4/3 * e * sqrt(r) * pow(penDepth, 3./2.) / dr : 0.;
          
          const SoAFloatPrecision fx = drx * fac;
          const SoAFloatPrecision fy = dry * fac;
          const SoAFloatPrecision fz = drz * fac;

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;
          if (newton3) {
            fxArr[j] = fx;
            fyArr[j] = fy;
            fzArr[j] = fz;
          }
        }
        
        // scatter the forces to where they belong, this is only needed for newton3
        if (newton3) {
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            const size_t j = neighborListPtr[joff + tmpj];
            fxptr[j] -= fxArr[tmpj];
            fyptr[j] -= fyArr[tmpj];
            fzptr[j] -= fzArr[tmpj];
          }
        }
      }
    }
    // this loop goes over the remainder and uses no optimizations
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;

      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision drx = xptr[indexFirst] - xptr[j];
      const SoAFloatPrecision dry = yptr[indexFirst] - yptr[j];
      const SoAFloatPrecision drz = zptr[indexFirst] - zptr[j];

      const SoAFloatPrecision drx2 = drx * drx;
      const SoAFloatPrecision dry2 = dry * dry;
      const SoAFloatPrecision drz2 = drz * drz;

      const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
      const SoAFloatPrecision dr = sqrt(dr2);

      const SoAFloatPrecision radI = radptr[indexFirst];
      const SoAFloatPrecision radJ = radptr[j];
      const SoAFloatPrecision rad = radI + radJ;

      if(dr >= rad | dr == 0) {continue;}

      const SoAFloatPrecision poissonI = poissonptr[indexFirst];
      const SoAFloatPrecision youngI = youngptr[indexFirst];
      const SoAFloatPrecision poissonJ = poissonptr[j];
      const SoAFloatPrecision youngJ = youngptr[j];

      const SoAFloatPrecision poissonI2 = poissonI * poissonI;
      const SoAFloatPrecision poissonJ2 = poissonJ * poissonJ;
      const SoAFloatPrecision e1 = (1. - poissonI2)/youngI;
      const SoAFloatPrecision e2 = (1. - poissonJ2)/youngJ;
      const SoAFloatPrecision esum = e1 + e2;
      const SoAFloatPrecision e = 1. / esum;

      const SoAFloatPrecision penDepth = rad - dr;

      const SoAFloatPrecision r = radI * radJ / rad;
      const SoAFloatPrecision fac = 4/3 * e * sqrt(r) * pow(penDepth, 3./2.) / dr;

      const SoAFloatPrecision fx = drx * fac;
      const SoAFloatPrecision fy = dry * fac;
      const SoAFloatPrecision fz = drz * fac;

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;
      
      if (newton3) {
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;
      }

    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;
    }
  }

    // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

};
}  // namespace autopas
