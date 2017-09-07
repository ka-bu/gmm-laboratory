#ifndef DESCRIPTION_H
#define DESCRIPTION_H

#include "../base.h"


class Description
{
	
public: 
	
	std::string getNameTag() const;

 	static void parseDescriptionString(Description* description, std::string descriptionstring);
	
	/**
	 * converts strings of the form a-b.c to three unsigned ints a, b and c, where b and c are optional
	 */
	static void readRangeSpec(std::string const& spec, unsigned int& from, unsigned int& to, unsigned int& step);
	
	/**
	 * converts strings of the form <range spec>,<range spec>,... to a std::vector<unsigned int>
	 */
	static std::vector<unsigned int> readIntList(std::istringstream& iss);
	
	/**
	 * converts strings of the form fp_type,fp_type, ... to a std::vector<fp_type>
	 */
	static std::vector<fp_type> readFpTypeList(std::istringstream& iss);
	
	static fp_type readFpType(std::istringstream& iss);
	
	static std::string readString(std::istringstream& iss);
	
	static std::string readUpperString(std::istringstream& iss);
	
	static unsigned int readInt(std::istringstream& iss);
	
	
protected:

	std::string nametag = "NoName";
	
};

#endif