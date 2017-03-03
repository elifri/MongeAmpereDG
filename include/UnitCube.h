/*
 * UnitCube.hh
 *
 *  Created on: Mar 5, 2015
 *      Author: friebel
 */

#ifndef SRC_UNITCUBE_HH_
#define SRC_UNITCUBE_HH_

# include <array>
# include <memory>

#include <config.h>
#include <vector>
# include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#ifdef HAVE_ALUGRID
# include <dune/grid/alugrid.hh>
# elif HAVE_DUNE_ALUGRID
# include <dune/alugrid/grid.hh>
# endif

using namespace Dune;

template <typename T>
struct Class {
    Class () {
        std::cout << "Default constructor" << std::endl;
    }
    template <typename U>
    Class (const Class<U>& rhs) {
        std::cout << "Const copy constructor" << std::endl;
    }
    template <typename U>
    Class (Class<U>& rhs)
        : Class (const_cast<const Class<U>&> (rhs))
    {
        std::cout << "Copy constructor (templated)" << std::endl;
    }
/* // DOES NOT WORK WITHOUT THE NEXT:
    Class (Class& rhs)
        : Class (const_cast<const Class&> (rhs))
    {
        std::cout << "Copy constructor (not templated)" << std::endl;
    }
*/
};



// default implementation for any template parameter
template<typename T>
class UnitCube {
public:
	typedef T GridType;

	static const int dim = GridType::dimension;

	// constructor throwing exception
	UnitCube() {
		Dune::FieldVector<typename GridType::ctype, dim> lowerLeft(0);
		Dune::FieldVector<typename GridType::ctype, dim> upperRight(1);
		std::array<unsigned int, dim> elements;
		std::fill(elements.begin(), elements.end(), 1);
		grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(
				lowerLeft, upperRight, elements);
	}

	T & grid() {
		return *grid_;
	}

private:
	//the constructed grid object
	std::shared_ptr<T> grid_;
};

//AluGrid specialisation
template<int dim>
class UnitCube<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> > {
public:
	typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> GridType;
	typedef Dune::FieldVector<typename GridType::ctype, dim> SpaceType;
private:
	std::shared_ptr<GridType> grid_;

public:
	UnitCube(SpaceType lowerLeft, SpaceType upperRight, int n) {
		std::array<unsigned int, dim> elements;
		std::fill(elements.begin(), elements.end(), 2);

//		grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(
//				lowerLeft, upperRight, elements);

		GridFactory<GridType> factory;

	    if(MPIHelper::getCollectiveCommunication().rank() == 0)
	    {
	        // Insert uniformly spaced vertices
	        std::array<unsigned int,dim> vertices = elements;
	        for (std::size_t i=0; i<vertices.size(); i++)
	          vertices[i]++;

	        // Create vertices
			FieldVector<double, GridType::dimensionworld> pos(0);
			pos[0] = lowerLeft[0];
			pos[1] = lowerLeft[1];
			factory.insertVertex(pos);

			pos[0] = upperRight[0];
			pos[1] = lowerLeft[1];
			factory.insertVertex(pos);

			pos[0] = (lowerLeft[0]+upperRight[0])/2.;
			pos[1] = (lowerLeft[1]+upperRight[1])/2.;
			factory.insertVertex(pos);


			pos[0] = lowerLeft[0];
			pos[1] = upperRight[1];
			factory.insertVertex(pos);

			pos[0] = upperRight[0];
			pos[1] = upperRight[1];
			factory.insertVertex(pos);


	        // Insert the elements
	        std::vector<unsigned int> corners(dim+1);

	        corners[0] = 0;
	        corners[1] = 2;
	        corners[2] = 1;

            factory.insertElement
	              (GeometryType(GeometryType::simplex, dim),
	              corners);

	        // 'base' is the index of the lower left element corner
	        corners[0] = 1;
	        corners[1] = 2;
	        corners[2] = 4;

            factory.insertElement
	              (GeometryType(GeometryType::simplex, dim),
	              corners);

	        // 'base' is the index of the lower left element corner
	        corners[0] = 4;
	        corners[1] = 2;
	        corners[2] = 3;

            factory.insertElement
	              (GeometryType(GeometryType::simplex, dim),
	              corners);

	        corners[0] = 3;
	        corners[1] = 2;
	        corners[2] = 0;

            factory.insertElement
	              (GeometryType(GeometryType::simplex, dim),
	              corners);
	      }
	    grid_ = shared_ptr<GridType>(factory.createGrid());

		grid_->globalRefine(n);
	}

	enum crossingType {LeftDiagonal, RightDiagonal};

  UnitCube(SpaceType lowerLeft, SpaceType upperRight, int n, crossingType crossing) {
    std::array<unsigned int, dim> elements;
    std::fill(elements.begin(), elements.end(), 2);

//    grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(
//        lowerLeft, upperRight, elements);

    GridFactory<GridType> factory;

      if(MPIHelper::getCollectiveCommunication().rank() == 0)
      {
          // Insert uniformly spaced vertices
          std::array<unsigned int,dim> vertices = elements;
          for (std::size_t i=0; i<vertices.size(); i++)
            vertices[i]++;

              // Create vertices
          FieldVector<double, GridType::dimensionworld> pos(0);
          pos[0] = lowerLeft[0];
          pos[1] = lowerLeft[1];
          factory.insertVertex(pos);

          pos[0] = upperRight[0];
          pos[1] = lowerLeft[1];
          factory.insertVertex(pos);

          pos[0] = lowerLeft[0];
          pos[1] = upperRight[1];
          factory.insertVertex(pos);

          pos[0] = upperRight[0];
          pos[1] = upperRight[1];
          factory.insertVertex(pos);

          // Insert the elements
          std::vector<unsigned int> corners(dim+1);
          switch(crossing)
          {
          case LeftDiagonal:

            corners[0] = 0;
            corners[1] = 1;
            corners[2] = 3;

            factory.insertElement
                (GeometryType(GeometryType::simplex, dim),
                corners);

            // 'base' is the index of the lower left element corner
            corners[0] = 3;
            corners[1] = 2;
            corners[2] = 0;

            factory.insertElement
                (GeometryType(GeometryType::simplex, dim),
                corners);
            break;
          case RightDiagonal:
            corners[0] = 0;
            corners[1] = 1;
            corners[2] = 2;

            factory.insertElement
                (GeometryType(GeometryType::simplex, dim),
                corners);

            // 'base' is the index of the lower left element corner
            corners[0] = 3;
            corners[1] = 2;
            corners[2] = 1;

            factory.insertElement
                (GeometryType(GeometryType::simplex, dim),
                corners);
            break;
          }
        }
      grid_ = shared_ptr<GridType>(factory.createGrid());

    grid_->globalRefine(n);
  }


	UnitCube(int n) : UnitCube(SpaceType(0), SpaceType(1), n) {}

	GridType& grid() {
		return *grid_;
	}

	const std::shared_ptr<GridType>& grid_ptr() const{
		return grid_;
	}

};

template<int dim>
class UnitCube<Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> >> {
public:
	typedef Dune::YaspGrid<dim, EquidistantOffsetCoordinates<double,dim> > GridType;
  typedef Dune::FieldVector<typename GridType::ctype, dim> SpaceType;

  UnitCube(int size) {
		Dune::FieldVector<double, dim> length(1.0);
		std::array<int, dim> elements;
		std::fill(elements.begin(), elements.end(), size);

		grid_ = std::unique_ptr < GridType
				> (new GridType(length, elements));
	}

  UnitCube(SpaceType lowerLeft, SpaceType upperRight, int n)
  {
    std::array<int, dim> elements;
    std::fill(elements.begin(), elements.end(), n);

    grid_ = std::unique_ptr < GridType
        > (new GridType(lowerLeft, upperRight, elements));

  }

  const std::shared_ptr<GridType>& grid_ptr() const{
    return grid_;
  }

	GridType& grid() {
		return *grid_;
	}

private:
	//the constructed grid object
	std::shared_ptr<GridType> grid_;
};

//
#endif /* SRC_UNITCUBE_HH_ */
