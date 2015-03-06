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

#include <dune/grid/yaspgrid.hh>

#include <dune/alugrid/grid.hh>

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

private:
	std::shared_ptr<GridType> grid_;

public:
	UnitCube(int n) {
		Dune::FieldVector<typename GridType::ctype, dim> lowerLeft(0);
		Dune::FieldVector<typename GridType::ctype, dim> upperRight(1);
		std::array<unsigned int, dim> elements;
		std::fill(elements.begin(), elements.end(), n);

		grid_ = Dune::StructuredGridFactory<GridType>::createSimplexGrid(
				lowerLeft, upperRight, elements);
	}

	GridType& grid() {
		return *grid_;
	}
};

template<int dim>
class UnitCube<Dune::YaspGrid<dim> > {
public:
	typedef Dune::YaspGrid<dim> GridType;

	UnitCube(int size) {
		Dune::FieldVector<double, dim> length(1.0);
		std::array<int, dim> elements;
		std::fill(elements.begin(), elements.end(), size);

		grid_ = std::unique_ptr < Dune::YaspGrid<dim>
				> (new Dune::YaspGrid<dim>(length, elements));
	}

	GridType& grid() {
		return *grid;
	}

private:
	//the constructed grid object
	std::unique_ptr<GridType> grid_;
};

//
#endif /* SRC_UNITCUBE_HH_ */
