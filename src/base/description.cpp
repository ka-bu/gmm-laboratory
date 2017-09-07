#include "description.h"

#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>


std::string Description::getNameTag() const
{
	return this->nametag;
}

void Description::parseDescriptionString(Description* descriptions, std::string descriptionstring)
{
	std::cout << "Don't know how to parse " << descriptionstring << " yet." << std::endl;
	gmmlab_throw("Not yet implemented.");
}


void Description::readRangeSpec(std::string const& spec, unsigned int& from, unsigned int& to, unsigned int& step)
{
	std::istringstream subiss(spec);
	std::string entry;
	
	getline(subiss, entry, '-');
	from = atoi(entry.c_str());
	entry.clear();
	
	getline(subiss, entry, '.');
	if (entry.empty())
		to = from;
	else
		to = atoi(entry.c_str());
	entry.clear();
	
	getline(subiss, entry);
	if (entry.empty())
		step = 1;
	else
		step = atoi(entry.c_str());
}


std::vector<unsigned int> Description::readIntList(std::istringstream& iss)
{
	std::string list;
	getline(iss, list, ':');
	
	std::istringstream subiss(list);
	std::string entry;
	std::vector<unsigned int> v;	
	
	while (getline(subiss, entry, ','))
	{
		unsigned int from, to, step;
		readRangeSpec(entry, from, to, step);
		for (unsigned int i=from; i<=to; i+=step)
			v.push_back(i);
	}
	
	if(v.size()==0)
		gmmlab_throw("Empty List.");
	
	return v;
}

std::vector<fp_type> Description::readFpTypeList(std::istringstream& iss)
{
	std::string list;
	getline(iss, list, ':');
	
	std::istringstream subiss(list);
	std::string entry;
	std::vector<fp_type> v;
	
	while (getline(subiss, entry, ','))
		v.push_back(atof(entry.c_str()));
	
	if(v.size()==0)
		gmmlab_throw("Empty List.");
	
	return v;
}

fp_type Description::readFpType(std::istringstream& iss)
{
	std::string entry;
	getline(iss, entry, ':');
	if(entry.size()==0)
		gmmlab_throw("No FpType.")
	return atof(entry.c_str());
}

unsigned int Description::readInt(std::istringstream& iss)
{
	std::string entry;
	getline(iss,entry, ':');
	if(entry.size()==0)
		gmmlab_throw("No Int.")
	return atoi(entry.c_str());
}

std::string Description::readUpperString(std::istringstream& iss)
{
	std::string entry;
	getline(iss, entry, ':');
	boost::to_upper(entry);
	return entry;
}

std::string Description::readString(std::istringstream& iss)
{
	std::string entry;
	getline(iss, entry, ':');
	return entry;
}