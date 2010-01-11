#ifndef UTILS_HPP
#define UTILS_HPP

#define MDINLINE inline

/** Calculates the SQuaRe of 'double' x, returning 'double'. */

MDINLINE double SQR(double x) { return x*x; }

/** Returns the distance between two positions squared and stores the
    distance vector pos1-pos2 in vec.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 *  \param vec  vecotr pos1-pos2.
 *  \return distance squared
*/

MDINLINE double distance2vec(double pos1[3], double pos2[3], double vec[3])
{
  vec[0] = pos1[0]-pos2[0];
  vec[1] = pos1[1]-pos2[1];
  vec[2] = pos1[2]-pos2[2];
  return SQR(vec[0]) + SQR(vec[1]) + SQR(vec[2]);
}

#endif
