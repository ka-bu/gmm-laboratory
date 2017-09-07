#include "ioutil.h"

#include "../settings/settings.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdint>
#include <dirent.h>

const char IOUtil::SEPARATOR = ',';
const char IOUtil::COMMENT = '#';
const std::string IOUtil::LIST_FILE_EXTENSION = ".list";
const std::string IOUtil::NA = "";

IOUtil::IOUtil()
{
}

IOUtil::IOUtil(std::string file_extension)
{
	this->file_extension = file_extension;
}

IOUtil::~IOUtil()
{
	if(this->is_open())
		this->close();
}

void IOUtil::header()
{
}


void IOUtil::open(std::string const& path, unsigned int k)
{
	this->k = k;
	std::string name = path;
	name.append(this->file_extension);

	if(commonSettings().verbose)
		std::cout << "IOUtil::open() -- writing "<< this->file_extension <<" statistics to file " << name << std::endl;

	this->filestream.open(name.c_str());
	IOUtil::set_precision(this->filestream);
	if(!filestream.is_open())
	{
		std::ofstream tmp;
		tmp.open(name.c_str());
		tmp.close();
		this->filestream.open(name.c_str());
	}
	this->header();
	if(this->filestream.tellg()>0) // if a header has been added
		this->filestream << std::endl;
	
	if (!this->filestream.is_open())
	{
		std::cout << "IOUtil::open() -- File creation for path " << path << " failed!" << std::endl;
		gmmlab_throw("An IOUtil cannot be opened.")
	}
	
}

bool IOUtil::is_open()
{
	return this->filestream.is_open();
}

void IOUtil::close()
{
	if(this->filestream.is_open())
		this->filestream.close();
}

bool IOUtil::ends_with(std::string s, std::string postfix)
{
	return s.size()>=postfix.size() && s.substr(s.size()-postfix.size(), std::string::npos) == postfix;
}

void IOUtil::set_precision(std::fstream& filestream){
	filestream << std::setiosflags(std::ios::scientific)
		<< std::setiosflags(std::ios::showpoint)
		<< std::setprecision(std::numeric_limits<fp_type>::digits10);
}

void IOUtil::set_precision(std::ofstream& filestream){
	filestream << std::setiosflags(std::ios::scientific)
		<< std::setiosflags(std::ios::showpoint)
		<< std::setprecision(std::numeric_limits<fp_type>::digits10);
}


void IOUtil::getFilenames(std::string const& pattern, std::vector<std::string>& list)
{	
	// divide pattern into path and prefix
	std::string dirname, filepat;
	std::size_t index = pattern.find_last_of('/');
	if (index==std::string::npos)
		dirname = "./";
	else
		dirname = pattern.substr(0,index+1);
	filepat = pattern.substr(index+1,std::string::npos);
	
	// open directory
	if (commonSettings().verbose)
		std::cout << "searching directory " << dirname << std::endl;
	DIR* dir;
	dirent* entry;
	if ( (dir = opendir(dirname.c_str())) != NULL)
	{
		while ( (entry = readdir(dir)) != NULL)
			if (entry->d_name[0] != '.')
			{
				std::string filename(entry->d_name);
				if (!boost::filesystem::is_directory(boost::filesystem::path(dirname+"/"+filename)) &&  filename.substr(0,filepat.size())==filepat)
				{
					if (commonSettings().verbose)
						std::cout << "   found " << filename << std::endl;
					std::string path(dirname);
					path += filename;
					list.push_back(path);
				}
			}
		closedir(dir);
		if (commonSettings().verbose)
			std::cout << list.size() << " files found." << std::endl;
	}
	else
	{
		std::cout << "ioutil::getFilenames() - Trying to open directory " << dirname << std::endl;
		gmmlab_throw("Failed to open directory.");
	}
}


const std::string IOUtil::commonFileExtension(std::vector<std::string> const& list)
{
	if (list.empty())
		return "";

	std::size_t index = list.begin()->find_last_of('.');
	if (index==std::string::npos)
		return "";

	std::string extension = list.begin()->substr(index+1,std::string::npos);
	for (std::vector<std::string>::const_iterator it = list.begin(); it!=list.end(); ++it)
	{
		index = it->find_last_of('.');
		if (index==std::string::npos || it->substr(index+1,std::string::npos)!=extension)
			return "";
	}
	return extension;
}


std::vector<std::string> IOUtil::readDirectory( const std::string& path = std::string() )
{
	std::vector <std::string> result;

	// if a .list file is given, then read the filenames from there
	if(ends_with(path, LIST_FILE_EXTENSION))
	{
		std::ifstream dirfilestream;
		dirfilestream.open(path);
		std::string line;
		while(getline(dirfilestream, line))
		{
			if(line!="")
				result.push_back(line);
		}
	}
	// otherwise open the path, read the filenames and sort them in alphabetical order
	else
	{
		dirent* de;
		DIR* dp;
		errno = 0;
		dp = opendir(path.c_str());
		std::cout << "Directory: " << path << std::endl;

		if(!dp)
		{
			std::cout << "ioutil::readDirectory - Trying to open directory " << path << std::endl;
			gmmlab_throw("Cannot open directory");
		}
		else
		{
			while (true)
			{
				errno = 0;
				de = readdir( dp );
			  if (de == NULL) break;
			  result.push_back( std::string( de->d_name ) );
			}
			closedir( dp );
			std::sort( result.begin(), result.end() );
		}
	}
	return result;
}

void IOUtil::store(Vector const& entry)
{
	idx_type d = entry.size();
	for(idx_type i=0; i<d; ++i){
		filestream << entry(i);
		if(i<d-1)
			filestream << IOUtil::SEPARATOR;
	}
	if(d>0)
		filestream << std::endl;
}