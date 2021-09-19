/**
 * @file DEMFunctorPolyhedra.h
 *
 * @date 01.07.21
 * @author R. Penz
 */

#pragma once

#include <array>

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
 */

template <class Particle, FunctorN3Modes useNewton3 = FunctorN3Modes::Both>
class DEMFunctorPolyhedra
    : public Functor<Particle,
                     DEMFunctorPolyhedra<Particle, poisson1, poisson2, youngmod1, youngmod2, r1, r2>> {
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
  DEMFunctorPolyhedra() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit DEMFunctorPolyhedra(double cutoff, void * /*dummy*/)
      : Functor<Particle, DEMFunctorPolyhedra<Particle, poisson1, poisson2, youngmod1, youngmod2, r1, r2>>(),
        

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    if (i.isDummy() or j.isDummy()) {
      return;
    }

    //todo: AOS functor code DEMContact

    /**
     * 
     * equilateral Tetrahedron
     * 
     * needed Geometrical Attributes:
     * RotX, RotY, RotZ
     * CenterX, CenterY, CenterZ
     * Radius
     * Torque
     * 
     * non-equilateral Tetrahedron
     * 
     * * needed Geometrical Attributes:
     * CenterX
     * CenterY
     * CenterZ
     * V1x, V1y, V1z
     * V2x, V2y, V2z
     * V3x, V3y, V3z
     * V4x, V4y, V4z
     * Torque
     * 
     * 
     * Polyhedron equilateral/non-equilateral
     * 
     * multiple Tetrahedron together, complex Force Calculation required
     * 
     * Physical Attributes:
     * 
     * 
     */


    /**
     * Todo: imlement after tetrahedra functor/particle
     * 
     */

    auto dr = utils::ArrayMath::sub(i.getR(), j.getR()) //distance between ParticleCenters
    double penDepth = i.getRad()+j.getRad()-dr;

    if(penDepth <= 0) return; //return if Particles dont intersect

    double e1 = (1-pow(i.getPoisson(), 2))/i.getYoung()
    double e2 = (1-pow(j.getPoisson(), 2))/j.getYoung()
    double e = 1/(e1+e2);

    double r = i.getRad()*j.getRad()/(r1+r2);

    //calculate Force and ForceVector
    double f = 4/3*e*sqrt(r)*pow(penDepth, 3./2.);
    auto vecf = utils::ArrayMath::mulScalar(utils::ArrayMath::normalize(dr),f) 
    i.add(vecf);
    if(newton3)
    {
      j.sub(vecf);
    }
  }

  void SoAFunctorCalc(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2,  size_t shift)
  {
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles == 0) return;

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    const auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    const auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    const auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    const auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    const auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict rad1ptr = soa1.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson1ptr = soa1.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young1ptr = soa1.template begin<Particle::AttributeNames::young>();
    const auto *const __restrict rad2ptr = soa2.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poisson2ptr = soa2.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict young2ptr = soa2.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownership1ptr = soa.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownership2ptr = soa.template begin<Particle::AttributeNames::ownershipState>();
    
    for(unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      const auto ownedStateI = ownershipptr[i];
      if(ownedStateI == OwnershipState::dummy) {return;}

      //accumulating Force directions
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

    
// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
// shift for SoAFunctorSingle
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
     for(unsigned int j = i+shift; j < soa.getNumParticles(); ++j) {
        const auto ownedStateJ = ownershipptr[j];

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

        fx2ptr[j] -=fx;
        fy2ptr[j] -=fy;
        fz2ptr[j] -=fz;
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

    SoAFunctorCalc(soa, soa, 1);

  }

    /**
    * old version of SoAFunctorSingle, keep for the moment
    */

    /*
    if (soa.getNumParticles() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    const auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    const auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::force>();

    const auto *const __restrict radptr = soa.template begin<Particle::AttributeNames::rad>();
    const auto *const __restrict poissonptr = soa.template begin<Particle::AttributeNames::poisson>();
    const auto *const __restrict youngptr = soa.template begin<Particle::AttributeNames::young>();

    const auto *const __restrict ownershipptr = soa.template begin<Particle::AttributeNames::ownershipState>();
    
    for(unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      const auto ownedStateI = ownershipptr[i];
      if(ownedStateI == OwnershipState::dummy) {return;}

      //accumulating Force directions
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

    
// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
     for(unsigned int j = i+1; j < soa.getNumParticles(); ++j) {
        const auto ownedStateJ = ownershipptr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
        const SoAFloatPrecision dr = sqrt(dr2);

        const SoAFloatPrecision radI = radptr[i];
        const SoAFloatPrecision radJ = radptr[j];
        const SoAFloatPrecision rad = radI + radJ;

        // Mask away if particles arent intersecting or if j is dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy;

        const SoAFloatPrecision poissonI = poissonptr[i];
        const SoAFloatPrecision youngI = youngptr[i];
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
        const SoAFloatPrecision fac = mask ? 4/3 * e * sqrt(r) * pow(penDepth, 3./2.) / dr : 0.;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        fxptr[j] -=fx;
        fyptr[j] -=fy;
        fzptr[j] -=fz;
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
  }
  */

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
    SoaFunctorCalc(soa1, soa2, 0);
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
    return std::array<typename Particle::AttributeNames, 11>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::rad,    Particle::AttributeNames::poisson, 
        Particle::AttributeNames::young,  Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 8>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::rad,    Particle::AttributeNames::poisson, 
        Particle::AttributeNames::young,  Particle::AttributeNames::ownershipState};
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

          // const size_t j = currentList[jNeighIndex];

          const auto ownedStateJ = ownedStateArr[j];

          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
          const SoAFloatPrecision dr = sqrt(dr2);

          const SoAFloatPrecision radI = radptr[i];
          const SoAFloatPrecision radJ = radptr[j];
          const SoAFloatPrecision rad = radI + radJ;

          // Mask away if particles arent intersecting or if j is dummy.
          const bool mask = dr <= rad and ownedStateJ != OwnershipState::dummy;

          const SoAFloatPrecision poissonI = poissontmp[i];
          const SoAFloatPrecision youngI = youngtmp[i];
          const SoAFloatPrecision poissonJ = poissonarr[j];
          const SoAFloatPrecision youngJ = youngarr[j];

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

      if(dr <= rad) {continue};

      const SoAFloatPrecision poissonI = poissontmp[indexFirst];
      const SoAFloatPrecision youngI = youngtmp[indexFirst];
      const SoAFloatPrecision poissonJ = poissonarr[j];
      const SoAFloatPrecision youngJ = youngarr[j];

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

    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;
    }
  }

  // variables for AoSFunctor
  double _poisson1, _poisson2, _youngmod1, _youngmod2, _rad1, _rad2 = 0.;

};
}  // namespace autopas
