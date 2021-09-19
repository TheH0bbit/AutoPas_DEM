#include "DEMTest.h"

#include "autopas/AutoPas.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

#include <iostream>


TEST_F(DEMTest, GridFillwithBoxMin) {
  auto autoPas = autopas::AutoPas<DEMParticle>(std::cout);
  std::array<double, 3> boxmin = {5., 5., 5.};
  std::array<double, 3> boxmax = {10., 10., 10.};
  autoPas.setBoxMax(boxmax);
  autoPas.setBoxMin(boxmin);
  DEMParticle dummy;

  autoPas.init();
  autopasTools::generators::GridGenerator::fillWithParticles(autoPas, {5, 5, 5}, dummy, {1, 1, 1}, boxmin);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(autopas::utils::inBox(iter->getR(), boxmin, boxmax));
  }
}

//Test-Test
TEST_F(DEMTest, TestingTests) {
  EXPECT_TRUE(true);
}