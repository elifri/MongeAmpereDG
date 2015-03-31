/*
 * mixedElement.hpp
 *
 *  Created on: Mar 9, 2015
 *      Author: friebel
 */

#ifndef SRC_MIXEDELEMENT_HPP_
#define SRC_MIXEDELEMENT_HPP_

#include<memory>
#include<string>

struct u{
};
struct u_DH{
	std::string name;
};

template <class LocalFiniteElement0, class LocalFiniteElement1, class FirstVariable=u, class SecondVariable=u_DH >
class MixedElement{


private:
	enum ElementSwitch { LFUO=0, LFU1=1};

public:
	MixedElement()
	{
//		assert(false && "Error, you cannot initialise a mixed element without names");
//		exit(-1);
	}

	MixedElement(const std::string name0, const std::string name2)
	{
		name_to_element[name0] = LFUO;
		name_to_element[name0] = LFU1;
	}


/**
 * return the underlying geometry type for this elements
 * @return
 */
Dune::GeometryType type() const{
	assert (lfu0.type() == lfu1.type());
	return lfu0.type();
}

/**
 * return the number of dofs in this element
 * @return
 */
unsigned int size() const
{
	return lfu0.size()+lfu1.size();
}

unsigned int size(const FirstVariable &v0) const { return lfu0.size();}
unsigned int size(const SecondVariable &v0) const { return lfu1.size();}

unsigned int order() const
{
	return std::max(lfu0.localBasis().order(), lfu1.localBasis().order());
}

unsigned int order(const FirstVariable &v0) const { return lfu0.localBasis().order();}
unsigned int order(const SecondVariable &v0) const { return lfu1.localBasis().order();}


private:
	LocalFiniteElement0 lfu0;
	LocalFiniteElement1 lfu1;

public:
const LocalFiniteElement0* operator()(const FirstVariable &v0) const	{return &lfu0;}
const LocalFiniteElement1* operator()(const SecondVariable &v0) const {return &lfu1;}


private:
	std::map<std::string, ElementSwitch> name_to_element;

};




#endif /* SRC_MIXEDELEMENT_HPP_ */
