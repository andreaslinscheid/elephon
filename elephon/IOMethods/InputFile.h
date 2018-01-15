/*	This file InputFile.h is part of elephon.
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
 *  Created on: Nov 16, 2014
 *      Author: Andreas Linscheid
 */

#ifndef ELEPHON_IOMETHODS_INPUTFILE_H_
#define ELEPHON_IOMETHODS_INPUTFILE_H_

#include <string>
#include <vector>
#include <map>

namespace elephon
{
namespace IOMethods
{

/**
 * 	\brief Read input to the elephon code and parse it into a key/value map.
 *
 *	\ref parse_input performs that key value parsing of a string that contains the input.
 * 	Using the \ref get_input_config_value method one can check the input file for the value string that
 * 	corresponds to the key string. Internally, there is a mutable tracking of the keys accessed in this way.
 * 	One can \ref get_list_unread_input_parameters to examine at some point, which keys have
 * 	not been referenced so far. This allows to help the user to determine, say, misspelled input keys with
 * 	default values that are silently taken with their default otherwise.
 */
class InputFile
{
public:

	/**
	 * \brief Inspect the stored key/value pairs.
	 *
	 * Referenced keys are traced so that \ref get_list_unread_input_parameters() will not
	 * include the key "akey" any more if get_input_config_value( "akey" ) was called at some point.
	 *
	 * @param key String matching the key name.
	 * @return	The value string, or an empty string if there is no such key.
	 */
	std::string get_input_config_value( std::string const& key ) const;

	/**
	 * \brief Copy all the stored key/value pairs.
	 *
	 * Referenced keys are not traced so that \ref get_list_unread_input_parameters() will
	 * include the key "akey" any more if get_input_config_value( "akey" ) was called at some point.
	 *
	 * @return A map from the key strings to the value strings.
	 */
	std::map<std::string,std::string> get_all_input_config_values( ) const;

	/**
	 * \brief Get a list of keys in the input that have not been referenced.
	 *
	 * @return A vector with unreferenced keys.
	 */
	std::vector<std::string> get_list_unread_input_parameters() const;

	/**
	 * \brief Parse the content of the input file for 'key = value' string pairs.
	 *
	 *	See \ref parse_input on how the parsing of the content is done.
	 * @param fileName The name of the input file.
	 */
	void read_input_file(std::string const& fileName);

	/**
	 * \brief Parse the content of the input stream for 'key = value' string pairs.
	 *
	 *	See \ref parse_input on how the parsing of the content is done.
	 * @param s The input stream with the data.
	 */
	void read_input_stream(std::istream & s);

	/**
	 * \brief Parse the input string for 'key = value' string pairs.
	 *
	 * It ignores comments following the char '!' or '#' until the end of the line '\n'
	 * @param input The string where the key = value pairs are searched in.
	 */
	void parse_input(std::string const&input);
private:

	///Store the key/value pairs
	std::map<std::string,std::string> inputFileKeyValue_;

	///Store if get_input_config_value(key) was called at some point
	mutable std::map<std::string,bool> keyWasRead_;

	/**
	 * \brief Remove the string starting from '/', '!' or '#' until '\\n' from str.
	 * @param str Input string with not comment on output.
	 */
	void remove_comment(std::string & str) const;

	/**
	 * \brief Remove the leading and trailing whitespace chars from str.
	 * @param str Input string with not leading and trailing whitespace on output.
	 * @param whitespace The chars that count as white space.
	 */
	void trim_string(std::string & str, std::string const& whitespace = " \t\n") const;
};

} /* namespace IOMethods */
} /* namespace elephon */
#endif /* ELEPHON_IOMETHODS_INPUTFILE_H_ */
