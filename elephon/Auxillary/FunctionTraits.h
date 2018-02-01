/*	This file FunctionTraits.h is part of elephon.
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

#ifndef ELEPHON_AUXILLARY_FUNCTIONTRAITS_H_
#define ELEPHON_AUXILLARY_FUNCTIONTRAITS_H_

#include <type_traits>
#include <tuple>

namespace elephon
{
namespace Auxillary
{

/// \cond
namespace detail
{
template <typename T, bool TIsClass>
struct FunctionTraits_impl;
}
/// \endcond


/** Macro that generates the test for a given function name.
 *
 * 	generate a struct that perform the check if the first template paramter implements a
 * 	function of given name and (second template parameter) of given signature
 * 	using the SFINAE principle.
 */
#define MAKE_CHECK_HAS_FUNCTION_OF_THIS_NAME(name) 							                            \
																										\
template<typename, typename T>												                            \
struct has_##name {															                            \
	static_assert(															                            \
		std::integral_constant<T, false>::value,							                            \
		"Second template parameter needs to be of function type.");							            \
};																			                            \
                                                                                                        \
template<typename Obj, typename Ret, typename... Args>                                                  \
struct has_##name<Obj, Ret(Args...)> {                                            					    \
private:                                                                                                \
    template<typename T>                                                                                \
    static constexpr auto check(T*)                                                                     \
    	-> typename                                                                                     \
    	std::is_same<                                                                                   \
    			decltype(std::declval<T>().name( std::declval<Args>()... )),         					\
    			Ret                                                                                     \
    			>::type;                                                                                \
                                                                                                        \
    template<typename>                                                                                  \
    static constexpr std::false_type check(...);                                                        \
                                                                                                        \
    typedef decltype(check<Obj>(0)) type;                                                               \
                                                                                                        \
public:                                                                                                 \
                                                                                                        \
    static constexpr bool value = type::value;                                                          \
};

/**
 * Type trait to obtain the return type of a function or the operator() for classes and number of arguments.
 *
 * Usage: 	FunctionTraits<C>::nargs is an int with the number of arguemnts the opeartor() or function call accepts.
 * 			FunctionTraits<C>::return_type is the type of the return value.
 *
 * @tparam C	A class that implements the operator() or const variants or a function pointer.
 *
 */
template<class C>
struct FunctionTraits : private detail::FunctionTraits_impl<C,std::is_class<C>::value>
{
	/// The number of arguments the operator() or function accepts.
	static constexpr int nargs = detail::FunctionTraits_impl<C,std::is_class<C>::value>::nargs;

	/// The type of the return value of the operator() or the function.
	typedef typename detail::FunctionTraits_impl<C,std::is_class<C>::value>::result_type result_type;

	/// Obtain type information for argument i of the operator() or function.
	/// Usage as e.g. FunctionTraits<T>::template arg<0>::type
	template <std::size_t i>
	struct arg {
		typedef typename detail::FunctionTraits_impl<C,std::is_class<C>::value>::template arg<i>::type type;
	};
};

/// \cond
// Details of the implementation
namespace detail
{

//The following simply fills the body of every function trait object. Based on the specified templates
//	the number of arguments and their type, as well as the return type is recorded.
#define STRUCT_BODY_FUNCTION_TRAITS													\
static constexpr int nargs = sizeof...(Args);												\
																					\
typedef ResultType result_type;														\
																					\
template <std::size_t i> 																	\
struct arg																			\
{																					\
    typedef typename std::tuple_element<i, std::tuple<Args...> >::type type;		\
}

template <class T, bool TIsClass>
struct FunctionTraits_impl
{
	static_assert(
		std::integral_constant<T, false>::value,
		"First template parameter needs to be of Functor or function type.");
};

//The following specializations determine information about the operator() of a functor
template <class T>
struct FunctionTraits_impl<T,true> : public FunctionTraits_impl<decltype(&T::operator()),true> { };

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(Functor::*)(Args...) const,true> {
	typedef ResultType(Functor::*signature_type)(Args...) const;
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType const&(Functor::*)(Args...) const,true> {
	typedef ResultType const&(Functor::*signature_type)(Args...) const;
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType& (Functor::*)(Args...),true> {
	typedef ResultType& (Functor::*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(Functor::*)(Args...),true> {
	typedef ResultType(Functor::*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};


//The following specialization determine information about function pointers
template<typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(*)(Args...),false> {
	typedef ResultType (*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType& (*)(Args...),false> {
	typedef ResultType& (*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType* (*)(Args...),false> {
	typedef ResultType* (*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

//The following specialization determine information about member function pointers
template<class T,typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(T::*)(Args...),false> {
	typedef ResultType (T::*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class T,typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType& (T::*)(Args...),false> {
	typedef ResultType& (T::*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class T,typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType* (T::*)(Args...),false> {
	typedef ResultType* (T::*signature_type)(Args...);
	STRUCT_BODY_FUNCTION_TRAITS;
};

#undef STRUCT_BODY_FUNCTION_TRAITS
} /* namespace detail */
/// \endcond

} /* namespace Auxillary */
} /* namespace elephon */

#endif /* ELEPHON_AUXILLARY_FUNCTIONTRAITS_H_ */
