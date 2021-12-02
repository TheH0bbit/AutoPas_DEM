/**
 * @file TrapezGenerator.h
 * @author R. Penz
 * @date 11/09/21
 */

#pragma once

#include "autopas/utils/ParticleTypeTrait.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopasTools::generators {
/**
 * Generator for a one particle deep trapezoid of particles.
 */
class TrapezGenerator {
 public:
  /**
   * Fills any container (also AutoPas object) with a cuboid mesh of particles.
   * Particle properties will be used from the default particle. Particle IDs start from the default particle.
   * @tparam Container Arbitrary container class that needs to support addParticle().
   * @param container
   * @param particlesPerDim Number of particles per dimension.
   * @param defaultParticle
   * @param spacing Factor for distance between two particles along one dimension (default is 1).
   * @param offset Offset to move all particles.
   * @param inXDirection Direction in which the Trapez is built
   */
  template <class Container>
  static void fillWithParticles(Container &container, const std::array<size_t, 3> &particlesPerDim,
                                const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle =
                                    typename autopas::utils::ParticleTypeTrait<Container>::value(),
                                const bool &inXDirection = true,
                                const std::array<double, 3> &spacing = std::array<double, 3>{1., 1., 1.},
                                const std::array<double, 3> &offset = std::array<double, 3>{.5, .5, .5});
};

template <class Container>
void TrapezGenerator::fillWithParticles(
    Container &container, const std::array<size_t, 3> &dimensions,
    const typename autopas::utils::ParticleTypeTrait<Container>::value &defaultParticle,
    const bool &inXDirection,
    const std::array<double, 3> &spacing, const std::array<double, 3> &offset) {

size_t id = defaultParticle.getID();

if(inXDirection){ //not the most elegant solution, preferably change so that one can define the bottom and top row of the trapez and then have them connect

  //in x direction
  unsigned long length = dimensions[0];
  double baselength = (double)length * spacing[0];
  unsigned long depth = dimensions[1];
  unsigned long height = dimensions[2];
  double direction = 1;

  if(offset[0] < 0 and offset[1] < 0){direction = -1;}

//calculating the trapezs tilt/height to depth ratio
  double ratio = fabs((double)height/(double)depth);
  double yratio = sqrt((spacing[0]*spacing[0])/(ratio+1));
  int itcount = 0;
  size_t id = defaultParticle.getID();

  //main loop
  for (double y = yratio; y < depth; y = y + yratio)
    {
      for(double x = 0; x < baselength; x+=spacing[0]) 
      {
       //particles square
        auto p = defaultParticle;
        p.setR({x + offset[0], y*direction + offset[1], y*ratio + offset[2]});
        p.setID(id++);
       container.addParticle(p);
      }
      for(unsigned int a = 0; a < itcount; a++)
      {
        double x = (double)a * spacing[0];
        //particles left triangle
        auto p = defaultParticle;
        p.setR({offset[0] - x - spacing[0], y*direction + offset[1], y*ratio + offset[2]});
        p.setID(id++);
        container.addParticle(p);
        //particles right triangle
        p.setR({x + baselength + offset[0], y*direction + offset[1], y*ratio + offset[2]});
        p.setID(id++);
        container.addParticle(p);
      }
    itcount++;
    }
}
else{
  //in y direction
unsigned long length = dimensions[1];
double baselength = (double)length * spacing[0];
unsigned long depth = dimensions[0];
unsigned long height = dimensions[2];
double direction = 1;

if(offset[0] < 0 and offset[1] < 0){direction = -1;}

double ratio = fabs((double)height/(double)depth);
double yratio = sqrt((spacing[0]*spacing[0])/(ratio+1));
int itcount = 0;
  for (double x = yratio; x < depth; x+=yratio)
  {
    for(double y = 0; y < baselength; y+=spacing[0]) 
    {
      //particles square
      auto p = defaultParticle;
      p.setR({x*direction + offset[0], y + offset[1] , x*ratio + offset[2]});
      p.setID(id++);
      container.addParticle(p);
    }
    for(unsigned int a = 0; a < itcount; a++)
    {
      double y = (double)a * spacing[0];
      //particles left triangle
      auto p = defaultParticle;
      p.setR({x*direction + offset[0], offset[1] - y - spacing[0] , x*ratio + offset[2]});
      p.setID(id++);
      container.addParticle(p);
      //particles right triangle
      p.setR({x*direction + offset[0], y + baselength + offset[1], x*ratio + offset[2]});
      p.setID(id++);
      container.addParticle(p);
    }
    itcount++;
  }
}

}
}  // namespace autopasTools::generators
