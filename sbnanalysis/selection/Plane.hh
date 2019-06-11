#ifndef PLANE_H
#define PLANE_H

#include <string>
#include <iostream>
#include <vector>
#include "TVector3.h"

namespace selection{
  
  /**
   * @brief  Plane class
   */
  class Plane{

    public  :
  
      /**
       * @brief  Constructor
       *
       * @param  V, a point in the plane
       * @param  A, bounds of the plane
       * @param  B, bounds of the plane
       *
       */
      Plane(const TVector3 &V, const TVector3 &A, const TVector3 &B);

      // Getters
      /**
       * @brief  Get a point in the plane
       *
       * @return a point in the plane
       */
      TVector3 GetV() const;
      
      /**
       * @brief  Get a boundary
       *
       * @return one of the plane's boundaries
       */
      TVector3 GetA() const;

      /**
       * @brief  Get a boundary
       *
       * @return one of the plane's boundaries
       */
      TVector3 GetB() const;
      
      /**
       * @brief  Get the direction of a boundary
       *
       * @return one of the plane's boundary unit vectors
       */
      TVector3 GetUnitA() const;
      
      /**
       * @brief  Get the direction of a boundary
       *
       * @return one of the plane's boundary unit vectors
       */
      TVector3 GetUnitB() const;
      
      /**
       * @brief  Get the direction of the normal to the plane
       *
       * @return the normal to the plane
       */
      TVector3 GetUnitN() const;

      /**
       * @brief  Get the magnitude of the boundary width
       *
       * @return the magnitude of the boundary width
       */
      float    GetAlpha() const;
      
      /**
       * @brief  Get the magnitude of the boundary width
       *
       * @return the magnitude of the boundary width
       */
      float    GetBeta() const;

    private :
      
      // Member variables
      TVector3 m_a;     ///< Unit vector of the direction to a plane boundary
      TVector3 m_b;     ///< Unit vector of the direction to a plane boundary
      TVector3 m_n;     ///< Unit vector of the normal to the plane
      TVector3 m_A;     ///< Vector of the boundary width
      TVector3 m_B;     ///< Vector of the boundary width
      TVector3 m_V;     ///< Position within the plane
      float    m_alpha; ///< Magnitude of the boundary width
      float    m_beta;  ///< Magnitude of the boundary width
      
  }; // Plane
} // Selection
#endif
