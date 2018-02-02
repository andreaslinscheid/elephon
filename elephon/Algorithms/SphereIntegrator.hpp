/*	This file SphereIntegrator.hpp is part of elephon.
 *
 *  elephon is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  elephon is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with elephon.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Jan 31, 2018
 *      Author: A. Linscheid
 */

#include "Algorithms/SphereIntegrator.h"
#include "Algorithms/helperfunctions.hpp"

namespace elephon
{
namespace Algorithms
{

template<class Functor>
int
SphereIntegrator<Functor>::pick_rule_spherical_harmonic(int lMax) const
{
	std::vector<int> avail_rules = {3, 5, 7, 9, 11, 13, 15,
				17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77,
				83, 89, 95, 101, 107, 113, 119, 125, 131 };
	auto it = std::lower_bound(avail_rules.begin(), avail_rules.end(), lMax);
	return it == avail_rules.end() ? *avail_rules.rbegin() : *it;
}

template<class Functor>
int
SphereIntegrator<Functor>::get_num_pts_surface() const
{
	return numElem_;
}

template<class Functor>
void
SphereIntegrator<Functor>::get_surface_pts(
		elephon::Auxillary::alignedvector::aligned_vector<FAT> & thetas,
		elephon::Auxillary::alignedvector::aligned_vector<SAT> & phis) const
{
	thetas.assign(thetas_.begin(), thetas_.end());
	phis.assign(phis_.begin(), phis_.end());
};

template<class Functor>
void
SphereIntegrator<Functor>::initialize(int lebedev_rule)
{
	this->get_rule_data(lebedev_rule);
}

template<class Functor>
typename SphereIntegrator<Functor>::RT
SphereIntegrator<Functor>::integrate(
		Functor const & f,
		int lebedev_rule)
{
	this->initialize(lebedev_rule);

	Auxillary::alignedvector::aligned_vector<RT> data(numElem_);

	this->get_function_data(
			f,
			thetas_.data(),
			phis_.data(),
			static_cast<int>(thetas_.size()),
			data.data());

	assert(weights_.size() == data.size());
	auto itw = weights_.begin();
	RT integral = RT(0.0);
	for (auto it = data.begin(); it != data.end(); ++it, ++itw)
	{
		integral += (*it)*(*itw);
	}
	return integral*4.0*M_PI;
}

template<class Functor>
void
SphereIntegrator<Functor>::get_rule_data(int numRule)
{
	if (numRule_ == numRule)
		return;

	// update data
	numRule_ = numRule;

	// we attempt to format the chain of loadings below in a way that is most clear for the few numbers
	// that are shuffled around.
	std::vector<double> plainData;
	if ( 3 == numRule ) {
		numElem_ =  6; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_003.txt"
			};} else
	if ( 5 == numRule ) {
		numElem_ = 14; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_005.txt"
		};}	else
	if ( 7 == numRule ) {
		numElem_ = 26; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_007.txt"
			};} else
	if ( 9 == numRule) {
		numElem_ = 38; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_009.txt"
			};} else
	if ( 11 == numRule ){
		numElem_ = 50; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_011.txt"
			};} else
	if ( 13 == numRule ){
		numElem_ =74; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_013.txt"
			};} else
	if ( 15 == numRule ){
		numElem_ = 86; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_015.txt"
			};} else
	if ( 17 == numRule ){
		numElem_ = 110; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_017.txt"
			};} else
	if ( 19 == numRule ){
		numElem_ = 146; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_019.txt"
			};} else
	if ( 21 == numRule ){
		numElem_ = 170; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_021.txt"
			};} else
	if ( 23 == numRule ){
		numElem_ = 194; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_023.txt"
			};} else
	if ( 25 == numRule ){
		numElem_ =230; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_025.txt"
			};} else
	if ( 27 == numRule ) {
		numElem_ =  266; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_027.txt"
			};} else
	if ( 29 == numRule ) {
		numElem_ = 302; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_029.txt"
			};}	else
	if ( 31 == numRule ) {
		numElem_ = 350; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_031.txt"
			};} else
	if ( 35 == numRule) {
		numElem_ = 434; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_035.txt"
			};} else
	if ( 41 == numRule ){
		numElem_ = 590; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_041.txt"
			};} else
	if ( 47 == numRule ){
		numElem_ =770; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_047.txt"
			};} else
	if ( 53 == numRule ){
		numElem_ = 974; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_053.txt"
			};} else
	if ( 59 == numRule ){
		numElem_ = 1202; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_059.txt"
			};} else
	if ( 65 == numRule ){
		numElem_ = 1454; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_065.txt"
			};} else
	if ( 71 == numRule ){
		numElem_ = 1730; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_071.txt"
			};} else
	if ( 77 == numRule ){
		numElem_ = 2030; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_077.txt"
			};} else
	if ( 83 == numRule ){
		numElem_ =2354; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_083.txt"
			};} else
	if ( 89 == numRule ){
		numElem_ = 2702; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_089.txt"
			};} else
	if ( 95 == numRule ){
		numElem_ = 3074; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_095.txt"
			};} else
	if ( 101 == numRule ){
		numElem_ = 3074; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_101.txt"
			};} else
	if ( 107 == numRule ){
		numElem_ = 3890; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_107.txt"
			};} else
	if ( 113 == numRule ){
		numElem_ = 4334; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_113.txt"
			};} else
	if ( 119 == numRule ){
		numElem_ = 4802; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_119.txt"
			};} else
	if ( 125 == numRule ){
		numElem_ =5294; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_125.txt"
			};} else
	if ( 131 == numRule ){
		numElem_ =5810; plainData = std::vector<double>{
			#include "Algorithms/Lebedev_rule_data/lebedev_131.txt"
			};}
	else
	{
		throw std::logic_error("Invalid Lebedev rule specified");
	};

	if (numElem_*3 != static_cast<int>(plainData.size()))
	{
		throw std::logic_error("data loaded from file does not match the size expected in the program");
	}

	weights_.resize(numElem_);
	thetas_.resize(numElem_);
	phis_.resize(numElem_);

	// a little annoyance is that the data from the website uses different angle convention that what we do.
	// We simply take the points as defined on that website and convert the 3D cartesian coordinates into our convention.

	for (int id = 0 ; id < numElem_; ++id)
	{
		double thetaW = (plainData[id*3+0]/180.0)*M_PI;
		double phiW = (plainData[id*3+1]/180.0)*M_PI;
		weights_[id] = plainData[id*3+2];

		double x = std::cos(thetaW)*std::sin(phiW);
		double y = std::sin(thetaW)*std::sin(phiW);
		double z =                  std::cos(phiW);

		double r;
		Algorithms::helperfunctions::compute_spherical_coords(x,y,z, r,thetas_[id],phis_[id]);
		assert((r-1.0)<1e-8);
	}
}

namespace detail
{

MAKE_CHECK_HAS_FUNCTION_OF_THIS_NAME(evaluate_many);

template<class F, bool hasEvalMany>
struct eval {};

template<class F>
struct eval<F,true>
{
	// dispatch to evaluate many
	template<typename T1, typename T2, typename T3>
	static void evaluate( F const & f, T1 const * theta, T2 const * phi, int numEvals, T3 * data)
	{
		f.evaluate_many(theta, phi, numEvals, data);
	}
};

template<class F>
struct eval<F,false>
{
	template<typename T1, typename T2, typename T3>
	static void evaluate( F const & f, T1 const * theta, T2 const * phi, int numEvals, T3 * data)
	{
		for (int id = 0 ; id < numEvals; ++id)
			data[id] = f(theta[id], phi[id]);
	}
};
};

template<class Functor>
void
SphereIntegrator<Functor>::get_function_data(
		Functor const & f,
		FAT const * theta,
		SAT const * phi,
		int numEvals,
		RT * data) const
{
	detail::eval<Functor,
				detail::has_evaluate_many<Functor, void(FAT const*, SAT const*,int, RT*)>::value
				>
		::evaluate(f, theta, phi, numEvals, data);
}

} /* namespace Algorithms */
} /* namespace elephon */
