//  (C) Copyright Thomas Brandes, SCAI Fraunhofer

// Boost.Test

// Attention: #include <boost/test/unit_test.hpp> does not include main program 

#include <boost/test/included/unit_test.hpp>

using boost::unit_test_framework::test_suite;

#include <ParticleContainer.hpp>
#include <BasicProperty.hpp>

BOOST_AUTO_TEST_CASE( ParticleContainerTest ) {

    int bigBlocks = 10;
    int sizeBigBlocks = 15;

    ParticleContainer *container = new ParticleContainer(bigBlocks, sizeBigBlocks);

    ArrayProperty<int,3> position     ("position");

    container->addProperty (&position);

    int key = 1131;

    for (int block = 0; block < bigBlocks; block++) {

      for (int j = 0; j < sizeBigBlocks; j++) {

         ParticleRef refParticle = container->addParticle (block, key++);
 
         int (&pos) [3] = position [refParticle];

         pos [0] = 1;
         pos [1] = 2;
         pos [2] = -1;
     }
   }

   int sumpos [3];

   sumpos[0] = 0;
   sumpos[1] = 0;
   sumpos[2] = 0;

   for (ParticleIterator it = container->begin(); it !=container->end(); ++it) {

      ParticleRef refParticle = *it;

      sumpos[0] += position[refParticle][0];
      sumpos[1] += position[refParticle][1];
      sumpos[2] += position[refParticle][2];
   }

   BOOST_CHECK_EQUAL (sumpos[0], bigBlocks * sizeBigBlocks);
   BOOST_CHECK_EQUAL (sumpos[1], 2 * bigBlocks * sizeBigBlocks);
   BOOST_CHECK_EQUAL (sumpos[2], - bigBlocks * sizeBigBlocks);

}

test_suite*
init_unit_test_suite( int argc, char * argv[] ) {

    test_suite* test= BOOST_TEST_SUITE( "Particle Container Test" );

    return test;
}

