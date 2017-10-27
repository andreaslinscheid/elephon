/*	This file InputBase.h is part of elephon.
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
 *  Created on: Nov 15, 2014
 *      Author: Andreas Linscheid
 */

#ifndef ELEPHON_IOMETHODS_INPUTBASE_H_
#define ELEPHON_IOMETHODS_INPUTBASE_H_

#include "IOMethods/InputFile.h"
#include <string>
#include <vector>
#include <type_traits>
#include <string>
#include <stdexcept>

namespace elephon
{
namespace IOMethods
{

/**
 * 	\brief A CRTP class that is the base of all input classes in the code.
 *
 * 	Using the macros \ref INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT and
 * 	\ref INPUTBASE_INPUT_OPTION_MACRO together with this class it should be easy to
 * 	write a new input class. Just derive from this class as done in \ref input::test::InputBaseTest
 * 	and use a macro to add the input option.
 *
 * 	Consequently a manual entry is created and the method get_<VariableName>() can be used
 * 	to access the new option from the Input file.
 *
 */
template<class derived>
class InputBase
{

public:

	/**
	 * \brief default constructor.
	 */
	InputBase();

	/**
	 * \brief directly calls \ref parse_variables upon construction.
	 * @param inputFile
	 */
	InputBase(InputFile const& inputFile);

	/**
	 * \brief Write the internally stored input manual to the fileName.
	 *
	 * @param fileName The full path of the manual to be created.
	 */
	void build_input_manual(std::string const& fileName) const;

	/**
	 * \brief Take an input file and fill the internal variables with content.
	 *
	 * @param inputFile The parsed input file.
	 */
	void parse_variables(InputFile const& inputFile);

protected:

	///True if \ref parse_variables was called.
	bool isInit_ = false;

	/**
	 * \brief Set value to the transformed value string or defaultValue of input option type.
	 *
	 *	The defaultValue is used when valueString is empty.
	 *
	 * @param valueString The string that defines the value.
	 * @param defaultValue	The default if the value string is empty.
	 * @param value	The value which is defaultValue or the transformed valueString on output.
	 */
	template<typename T>
	void get_option(std::string const& valueString,T const& defaultValue,T &value) const;

	/**
	 * \brief Set value to the transformed non-empty value string of input option type.
	 *
	 * The method error out if valueString is empty.
	 *
	 * @param valueString The non-empty string that defines the value.
	 * @param value The transformed valueString on output.
	 */
	template<typename T>
	void get_option(std::string const& valueString,T &value) const;

	/**
	 * \brief Similar to \ref get_option() .
	 */
	template<typename T>
	void get_option(std::string const& valueString,std::vector<T> const& defaultValue,std::vector<T> &values) const;

	/**
	 * \brief Similar to \ref get_option() .
	 */
	template<typename T>
	void get_option(std::string const& valueString,std::vector<T> &values) const;

	/**
	 * \brief Add a manual entry with exactly this text.
	 *
	 * @param textWithOptionDescription The text to be appended to the manual.
	 */
	void add_option_to_manual(std::string const& textWithOptionDescription);

	/**
	 * Register an input option that will be parsed if \ref parse_variables() is called.
	 *
	 * Add the member function pointer of the Input class, derived from this, which parses
	 * some variable. The member function is of the type
	 *  void parse_dummyVarible(InputFile const&) that looks for "dummyVarible" in
	 *  the \ref InputFile and sets the value of _dummyVarible in the derived class
	 *  to the value that is found.
	 *
	 * @param function A member function ptr of the derived class that performs the parsing.
	 */
	void register_variable_parsing( void( derived:: * function )(InputFile const&) );

private:

	///Contains the full manual as set by the options in the derived class.
	std::string manual_;

	///Contains the  member function ptrs of the derived class that perform the parsing.
	std::vector< void(derived::*)(InputFile const&) > mebrFctnPtrToVariableParsing_;
};

/**
 * \brief Allows to add comma in macro expansion values.
 *
 * Use for example {0 COMMA_SUBSTITUTION 1 COMMA_SUBSTITUTION 2 COMMA_SUBSTITUTION 3}
 * to initialize, say, a vector with the list {0,1,2,3} via a macro variable
 */
#define COMMA_SUBSTITUTION ,

/**
 * \brief This macro allows to easily add new input variables with a default value.
 *
 *	The strategy is that all body of the derived class is constructed from this macro.
 *	If a new option is included, we need 1) store the variable's value internally and 2)
 *	have a method that return the value from the Input object. 3) We also need to
 *	register the parsing of the corresponding variable so that it will be read from the input.
 *	4) Furthermore would we like to generate a manual entry from the description.
 *	5) The directly derived class should be allowed to process, i.e. to change input variables.
 *	This require non-locality. Thus, an auxiliary void pointer is included that is initialized
 *	with a function call_##quantityName##_processing(). Inside this function (which returns the Nullptr)
 *	we add the manual entry and register the parsing of the variable.
 *	To add a new input variable use
 * 		INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(
 * 			quantityName,
 * 			description,
 * 			defaultValueText,
 * 			defaultValueStatement,
 * 			typeQ);
 *
 *	Here:
 *\param quantityName  		Is the name of quantity in the input file and
 *							internally in the Input module
 *\param description		A description of the variable that appears in the manual.
 *\param defaultValueText	A description text of the default value such as "zero" or "0" or "sequence 0 to 3".
 *\param defaultValueStatement	A compiling expression for the default such as 0 or
 *							{0 COMMA_SUBSTITUTION 1 COMMA_SUBSTITUTION 2 COMMA_SUBSTITUTION 3}
 *\param typeQ				The type of the variable such as double or size_t or std::vector<size_t>
 *
 */
#define INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT(											\
		quantityName,description,defaultValueText,defaultValueStatement,typeQ)				\
public:                                                                                 	\
	typeQ const& get_##quantityName() const {                                            	\
		if ( not this->isInit_ )															\
		{																					\
			throw std::logic_error("Calling get_"#quantityName"()"							\
									"before calling parse_variables()");					\
		};																					\
		return _##quantityName;                                                         	\
	};                                                                                  	\
public:                                                                                 	\
	void set_##quantityName(typeQ val) {          		                                  	\
		if ( not this->isInit_ )															\
		{																					\
			throw std::logic_error("Calling set_"#quantityName"()"							\
									"before calling parse_variables()");					\
		};																					\
		_##quantityName = val;                                                         		\
	};                                                                                  	\
private:																					\
	typeQ _##quantityName;																	\
	void parse_variable_##quantityName(::elephon::IOMethods::InputFile const& inputFile) {	\
		typeQ defaultValue = defaultValueStatement;											\
		std::string valueInInputFile = inputFile.get_input_config_value(#quantityName);		\
		this->get_option(valueInInputFile,defaultValue, _##quantityName);					\
	}																						\
	void * perform_##quantityName##_processing() {											\
		this->add_option_to_manual( std::string()+											\
		  "\n========================================================================\n"	\
			"||Variable :             || "+#quantityName+"\n"								\
			"||-----------------------||---------------------------------------------\n"	\
			"||Type is :              || "+#typeQ+"\n"										\
			"||-----------------------||---------------------------------------------\n"	\
			"||The default value is : || "+defaultValueText+"\n"							\
			"========================================================================\n"	\
			"|Description:|\n"																\
			"--------------\n"																\
			""+description+"\n"                                      						\
			"========================================================================\n\n");\
		this->register_variable_parsing(													\
			&std::remove_pointer<decltype(this)>::type ::parse_variable_##quantityName);	\
		return 0;																			\
	};																						\
	void * call_##quantityName##_processing													\
				= this->perform_##quantityName##_processing()								\

/**
 * \brief Similar macro as \ref INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT without default values.
 *
 * 	This requires the variable to be set in the input process or the program errors out.
 * 	See documentation of \ref INPUTBASE_INPUT_OPTION_MACRO_WITH_DEFAULT
 */
#define INPUTBASE_INPUT_OPTION_MACRO(quantityName,description,typeQ)						\
public:                                                                                 	\
	typeQ const& get_##quantityName() const {                                            	\
		if ( not this->_isInit ){															\
			throw std::runtime_error("Calling get_"#quantityName"()"						\
									"before calling parse_variables()",1);					\
		};																					\
		return _##quantityName;                                                         	\
	};                                                                                  	\
public:                                                                                 	\
	void set_##quantityName(typeQ val) {             		                               	\
		if ( not this->isInit_ )															\
		{																					\
			throw std::logic_error("Calling set_"#quantityName"()"							\
									"before calling parse_variables()");					\
		};																					\
		_##quantityName = val;                                                         		\
	};                                                                                  	\
private:																					\
	typeQ _##quantityName;																	\
	void parse_variable_##quantityName(::elephon::IOMethods::InputFile const& inputFile) {	\
		std::string valueInInputFile = inputFile.get_input_config_value(#quantityName);		\
		this->get_option(valueInInputFile,_##quantityName);									\
	}																						\
	void * perform_##quantityName##_processing() {											\
		this->add_option_to_manual( std::string()+											\
		  "\n========================================================================\n"	\
			"||Variable :             || "+#quantityName+"\n"								\
			"||-----------------------||---------------------------------------------\n"	\
			"||Type is :              || "+#typeQ+"\n"										\
			"||-----------------------||---------------------------------------------\n"	\
			"||The default value :    || This variable is MANDATORY!\n"						\
			"========================================================================\n"	\
			"|Description:|\n"																\
			"--------------\n"																\
			""+description+"\n"                                      						\
			"========================================================================\n\n");\
		this->register_variable_parsing(													\
			&std::remove_pointer<decltype(this)>::type ::parse_variable_##quantityName);	\
		return 0;																			\
	};																						\
	void * call_##quantityName##_processing													\
				= this->perform_##quantityName##_processing()								\


} /* namespace IOMethods */
} /* namespace elephon */
#include "IOMethods/InputBase.hpp"
#endif /* ELEPHON_IOMETHODS_INPUTBASE_H_ */
